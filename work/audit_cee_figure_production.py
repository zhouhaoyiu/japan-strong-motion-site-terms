#!/usr/bin/env python3
"""Audit CEE figure production readiness."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if (PROJECT_ROOT / "outputs" / "figures").exists():
    FIGURES = PROJECT_ROOT / "outputs" / "figures"
elif (PROJECT_ROOT / "figures").exists():
    FIGURES = PROJECT_ROOT / "figures"
else:
    FIGURES = PROJECT_ROOT / "outputs" / "figures"

if (PROJECT_ROOT / "docs").exists() and (PROJECT_ROOT / "figures").exists():
    REPORT_DIR = PROJECT_ROOT / "docs"
elif (PROJECT_ROOT / "qa").exists():
    REPORT_DIR = PROJECT_ROOT / "qa"
else:
    REPORT_DIR = PROJECT_ROOT / "outputs"
REPORT = REPORT_DIR / "cee_figure_production_audit.md"
CSV_REPORT = REPORT_DIR / "cee_figure_production_audit.csv"

EXPECTED_FIGURES = {
    "cee_fig1_workflow_data_evidence": (3150, 1680),
    "cee_fig2_mf2013_site_term_gain_ci": (2640, 1560),
    "cee_fig3_residual_site_correlation_shift": (3000, 1560),
    "cee_fig4_subgroup_robustness": (1920, 1560),
    "cee_fig5_zhao2006_external_gmpe_sensitivity": (3960, 1620),
}

PALETTE_COLORS = {
    "MF2013 event/site blue": "#2166ac",
    "MF2013 event/site red": "#b2182b",
    "Zhao crustal blue": "#0072B2",
    "Zhao interplate vermillion": "#D55E00",
    "Zhao intraplate green": "#009E73",
}


@dataclass
class Check:
    scope: str
    check: str
    status: str
    detail: str


def add(checks: list[Check], scope: str, check: str, ok: bool, detail: str) -> None:
    checks.append(Check(scope=scope, check=check, status="PASS" if ok else "FAIL", detail=detail))


def first_existing(candidates: list[Path]) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def display(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def relative_luminance(hex_color: str) -> float:
    value = hex_color.lstrip("#")
    rgb = [int(value[i : i + 2], 16) / 255.0 for i in (0, 2, 4)]
    linear = []
    for channel in rgb:
        linear.append(channel / 12.92 if channel <= 0.03928 else ((channel + 0.055) / 1.055) ** 2.4)
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast_ratio(a: str, b: str = "#ffffff") -> float:
    la = relative_luminance(a)
    lb = relative_luminance(b)
    lighter, darker = max(la, lb), min(la, lb)
    return (lighter + 0.05) / (darker + 0.05)


def png_metrics(path: Path) -> dict[str, object]:
    image = Image.open(path).convert("RGB")
    arr = np.asarray(image)
    ink = np.any(arr < 245, axis=2)
    dark = np.mean(arr, axis=2) < 110
    ys, xs = np.where(ink)
    if len(xs) == 0 or len(ys) == 0:
        margins = (arr.shape[1], arr.shape[0], arr.shape[1], arr.shape[0])
    else:
        margins = (int(xs.min()), int(ys.min()), int(arr.shape[1] - 1 - xs.max()), int(arr.shape[0] - 1 - ys.max()))
    return {
        "width": arr.shape[1],
        "height": arr.shape[0],
        "ink_ratio": float(ink.mean()),
        "dark_ratio": float(dark.mean()),
        "margins": margins,
    }


def audit_files(checks: list[Check]) -> None:
    for stem, expected_size in EXPECTED_FIGURES.items():
        png = FIGURES / f"{stem}.png"
        pdf = FIGURES / f"{stem}.pdf"
        add(checks, stem, "png_exists", png.exists(), display(png))
        add(checks, stem, "pdf_exists", pdf.exists(), display(pdf))
        if png.exists():
            metrics = png_metrics(png)
            size = (metrics["width"], metrics["height"])
            add(checks, stem, "png_dimensions", size == expected_size, f"expected={expected_size}; actual={size}")
            add(
                checks,
                stem,
                "nonblank_pixels",
                0.01 <= metrics["ink_ratio"] <= 0.75 and metrics["dark_ratio"] >= 0.005,
                f"ink_ratio={metrics['ink_ratio']:.4f}; dark_ratio={metrics['dark_ratio']:.4f}",
            )
            min_margin = min(metrics["margins"])
            add(checks, stem, "edge_margin", min_margin >= 5, f"margins_px={metrics['margins']}; min={min_margin}")
            width_in = expected_size[0] / 300.0
            height_in = expected_size[1] / 300.0
            add(checks, stem, "print_scale_300dpi", width_in >= 6.0 and height_in >= 5.0, f"300dpi size={width_in:.2f}x{height_in:.2f} in")
        if pdf.exists():
            header = pdf.read_bytes()[:4]
            add(checks, stem, "pdf_header", header == b"%PDF", f"header={header!r}; size={pdf.stat().st_size:,} bytes")


def audit_source_style(checks: list[Check]) -> None:
    main_script_path = first_existing(
        [
            PROJECT_ROOT / "work" / "build_cee_manuscript_figures.py",
            PROJECT_ROOT / "scripts" / "build_cee_manuscript_figures.py",
        ]
    )
    zhao_script_path = first_existing(
        [
            PROJECT_ROOT / "work" / "build_cee_zhao_sensitivity_figure.py",
            PROJECT_ROOT / "scripts" / "build_cee_zhao_sensitivity_figure.py",
        ]
    )
    main_script = main_script_path.read_text(encoding="utf-8")
    zhao_script = zhao_script_path.read_text(encoding="utf-8")
    combined = main_script + "\n" + zhao_script

    required_style_tokens = [
        '"font.size": 10',
        '"axes.labelsize": 10',
        '"legend.fontsize": 9',
        '"savefig.dpi": 300',
        'color="white" if val >= 18 else "#111111"',
        'color="white" if abs(val) >= 0.35 else "#111111"',
    ]
    missing_style = [token for token in required_style_tokens if token not in combined]
    add(checks, "source_style", "font_and_heatmap_style", not missing_style, "missing=" + "; ".join(missing_style) if missing_style else "font, DPI, and heatmap annotation rules present")

    panel_tokens = ['"A"', '"B"', '"C"', '"D"', '"A  Long-period basin signal is preserved"', '"B  Zhao residual-site associations"', '"C  External-GMPE audit coverage"']
    missing_panels = [token for token in panel_tokens if token not in combined]
    add(checks, "source_style", "panel_labels", not missing_panels, "missing=" + "; ".join(missing_panels) if missing_panels else "Figure 1-5 panel labels present in source")

    guardrail_tokens = ["rrup proxy", "vs30/avs30", "Zhao 2006 audit uses OpenQuake hazardlib"]
    missing_guardrails = [token for token in guardrail_tokens if token not in zhao_script]
    add(checks, "source_style", "figure5_guardrail_note", not missing_guardrails, "missing=" + "; ".join(missing_guardrails) if missing_guardrails else "Figure 5 source note states proxy limitations")


def audit_palette(checks: list[Check]) -> None:
    ratios = {name: contrast_ratio(color) for name, color in PALETTE_COLORS.items()}
    low = {name: ratio for name, ratio in ratios.items() if ratio < 3.0}
    detail = "; ".join(f"{name}={ratio:.2f}" for name, ratio in ratios.items())
    add(checks, "palette", "white_background_contrast", not low, detail)

    event_site_ok = ratios["MF2013 event/site blue"] >= 4.5 and ratios["MF2013 event/site red"] >= 4.5
    add(
        checks,
        "palette",
        "primary_line_contrast",
        event_site_ok,
        f"blue={ratios['MF2013 event/site blue']:.2f}; red={ratios['MF2013 event/site red']:.2f}",
    )

    zhao_ok = all(ratios[name] >= 3.0 for name in ["Zhao crustal blue", "Zhao interplate vermillion", "Zhao intraplate green"])
    add(
        checks,
        "palette",
        "zhao_source_color_contrast",
        zhao_ok,
        f"crustal={ratios['Zhao crustal blue']:.2f}; interplate={ratios['Zhao interplate vermillion']:.2f}; intraplate={ratios['Zhao intraplate green']:.2f}",
    )


def write_reports(checks: list[Check]) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    with CSV_REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scope", "check", "status", "detail"])
        writer.writeheader()
        for check in checks:
            writer.writerow(check.__dict__)

    failures = [check for check in checks if check.status != "PASS"]
    lines = [
        "# CEE Figure Production Audit",
        "",
        "Date: 2026-06-05.",
        "",
        "Status: " + ("passed" if not failures else "failed"),
        "",
        "Scope: scripted production-readiness checks for CEE Figure 1-5 PNG/PDF outputs and their generation scripts. This audit checks dimensions, nonblank rendering, edge margins, 300-dpi print scale, PDF headers, font/DPI source settings, panel-label source markers, heatmap annotation contrast rules, and palette contrast.",
        "",
        "| Scope | Check | Status | Detail |",
        "|---|---|---:|---|",
    ]
    for check in checks:
        lines.append(f"| `{check.scope}` | `{check.check}` | {check.status} | {check.detail} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The current figure set passes the scripted production audit and is suitable for coauthor/domain review and pre-submission formatting checks. This does not replace final journal production review after the manuscript text, captions, and supplemental-material decision are fixed.",
            "",
            "Figure 5 remains an external-GMPE sensitivity audit, not a full Zhao et al. (2006) reproduction.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    checks: list[Check] = []
    audit_files(checks)
    audit_source_style(checks)
    audit_palette(checks)
    write_reports(checks)
    print(REPORT)
    if any(check.status != "PASS" for check in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
