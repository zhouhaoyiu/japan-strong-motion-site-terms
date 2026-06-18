#!/usr/bin/env python3
"""Smoke-test the CEE v0.3 internal review package.

This checks package integrity and review evidence consistency without reading
raw J-SHIS/NIED flatfiles or local waveform/cache data.
"""

from __future__ import annotations

import hashlib
import py_compile
import shutil
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ZIP = PROJECT_ROOT / "outputs" / "cee_review_package_v0_3_internal.zip"
DEFAULT_EXTRACT_DIR = PROJECT_ROOT / "work" / "cee_review_package_v0_3_internal_smoke"
DEFAULT_REPORT = PROJECT_ROOT / "outputs" / "cee_review_package_v0_3_internal_smoke_test_report.md"

PACKAGE_DIRNAME = "cee_review_package_v0_3_internal"
FORBIDDEN_NAME_TOKENS = [
    "flatfile",
    "smrec_schema",
    "site_schema.tsv",
    "source_schema.tsv",
    ".mseed",
    ".BSON",
    "seisbench",
    "external_data",
    "__pycache__",
    ".pyc",
]
FORBIDDEN_TEXT_TOKENS = ["/Users/", "/private/", "/var/folders"]


@dataclass
class Check:
    name: str
    status: str
    detail: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def png_dimensions(path: Path) -> tuple[int, int]:
    header = path.read_bytes()[:24]
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path} is not a PNG")
    return struct.unpack(">II", header[16:24])


def add(checks: list[Check], name: str, ok: bool, detail: str) -> None:
    checks.append(Check(name=name, status="PASS" if ok else "FAIL", detail=detail))


def expected_file_checks(pkg: Path, checks: list[Check]) -> None:
    expected = [
        "README.md",
        "PACKAGE_MANIFEST.txt",
        "release/SHA256SUMS.txt",
        "release/environment.yml",
        "manuscript/cee_article_draft_v0_3.md",
        "manuscript/cee_combined_manuscript_v0_3_review.docx",
        "manuscript/cee_combined_manuscript_v0_3_review.pdf",
        "qa/docx_render_qa_report.md",
        "qa/package_integrity_qa_report.md",
        "tables/cee_claims_evidence_table.csv",
        "tables/jshis_mf2013_site_term_bootstrap_ci.csv",
        "tables/jshis_zhao2006_external_gmpe_station_site_correlations.csv",
        "tables/cee_fig5_zhao2006_external_gmpe_sensitivity_source.csv",
        "figures/cee_fig1_workflow_data_evidence.png",
        "figures/cee_fig2_mf2013_site_term_gain_ci.png",
        "figures/cee_fig3_residual_site_correlation_shift.png",
        "figures/cee_fig4_subgroup_robustness.png",
        "figures/cee_fig5_zhao2006_external_gmpe_sensitivity.png",
    ]
    missing = [item for item in expected if not (pkg / item).exists()]
    add(checks, "expected files", not missing, "missing: " + ", ".join(missing) if missing else f"{len(expected)} required files present")


def forbidden_scan(pkg: Path, checks: list[Check]) -> None:
    files = [path for path in pkg.rglob("*") if path.is_file()]
    name_hits = [
        str(path.relative_to(pkg))
        for path in files
        if any(token.lower() in str(path.relative_to(pkg)).lower() for token in FORBIDDEN_NAME_TOKENS)
    ]
    add(checks, "forbidden filenames", not name_hits, "hits: " + ", ".join(name_hits[:10]) if name_hits else "none")

    text_hits: list[str] = []
    for path in files:
        if path.suffix.lower() in {".md", ".py", ".txt", ".yml", ".yaml", ".csv"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if any(token in text for token in FORBIDDEN_TEXT_TOKENS):
                text_hits.append(str(path.relative_to(pkg)))
    add(checks, "local absolute paths", not text_hits, "hits: " + ", ".join(text_hits[:10]) if text_hits else "none")


def script_compile_checks(pkg: Path, checks: list[Check]) -> None:
    failures: list[str] = []
    scripts = sorted((pkg / "scripts").glob("*.py"))
    compile_dir = pkg.parent / "_py_compile"
    compile_dir.mkdir(parents=True, exist_ok=True)
    for script in scripts:
        try:
            py_compile.compile(str(script), cfile=str(compile_dir / f"{script.stem}.pyc"), doraise=True)
        except py_compile.PyCompileError as exc:
            failures.append(f"{script.name}: {exc.msg}")
    add(checks, "script compilation", not failures, "failures: " + "; ".join(failures) if failures else f"{len(scripts)} scripts compile")


def table_checks(pkg: Path, checks: list[Check]) -> None:
    tables = pkg / "tables"
    expectations = {
        "cee_claims_evidence_table.csv": (10, 7),
        "jshis_mf2013_site_term_bootstrap_ci.csv": (28, 19),
        "jshis_zhao2006_external_gmpe_station_site_correlations.csv": (32, 8),
        "cee_fig5_zhao2006_external_gmpe_sensitivity_source.csv": (31, 7),
        "jshis_zhao2006_external_gmpe_coverage_summary.csv": (4, 4),
    }
    failures: list[str] = []
    for name, expected_shape in expectations.items():
        df = pd.read_csv(tables / name)
        if df.shape != expected_shape:
            failures.append(f"{name}: expected {expected_shape}, got {df.shape}")
    add(checks, "table shapes", not failures, "failures: " + "; ".join(failures) if failures else f"{len(expectations)} expected table shapes match")

    boot = pd.read_csv(tables / "jshis_mf2013_site_term_bootstrap_ci.csv")
    all_event = boot[(boot["cluster_level"].eq("event")) & (boot["subgroup"].eq("all"))].copy()
    sa3 = all_event[all_event["target_label"].eq("SA(3.0s) RotD50")]["percent_mae_reduction"].iloc[0]
    sa1 = all_event[all_event["target_label"].eq("SA(1.0s) RotD50")]["percent_mae_reduction"].iloc[0]
    add(
        checks,
        "MF2013 key reductions",
        abs(sa1 - 16.067154) < 0.001 and abs(sa3 - 21.169213) < 0.001,
        f"SA1 event/all={sa1:.6f}; SA3 event/all={sa3:.6f}",
    )

    zhao = pd.read_csv(tables / "jshis_zhao2006_external_gmpe_station_site_correlations.csv")
    zhao_sa3_d1400 = zhao[(zhao["target_label"].eq("SA(3.0s) RotD50")) & (zhao["feature"].eq("log10(D1400)"))][
        "spearman_rho"
    ].iloc[0]
    add(
        checks,
        "Zhao D1400 key signal",
        abs(zhao_sa3_d1400 - 0.512581) < 0.001,
        f"SA3/D1400 rho={zhao_sa3_d1400:.6f}",
    )


def figure_and_document_checks(pkg: Path, checks: list[Check]) -> None:
    figures = pkg / "figures"
    expected_png_sizes = {
        "cee_fig1_workflow_data_evidence.png": (3150, 1680),
        "cee_fig2_mf2013_site_term_gain_ci.png": (2640, 1560),
        "cee_fig3_residual_site_correlation_shift.png": (3000, 1560),
        "cee_fig4_subgroup_robustness.png": (1920, 1560),
        "cee_fig5_zhao2006_external_gmpe_sensitivity.png": (3960, 1620),
    }
    failures: list[str] = []
    for name, expected_size in expected_png_sizes.items():
        size = png_dimensions(figures / name)
        if size != expected_size:
            failures.append(f"{name}: expected {expected_size}, got {size}")
        pdf = figures / name.replace(".png", ".pdf")
        if not pdf.read_bytes().startswith(b"%PDF"):
            failures.append(f"{pdf.name}: missing PDF header")
    add(checks, "figure files", not failures, "failures: " + "; ".join(failures) if failures else "5 PNG/PDF figure pairs verified")

    docx = pkg / "manuscript" / "cee_combined_manuscript_v0_3_review.docx"
    pdf = pkg / "manuscript" / "cee_combined_manuscript_v0_3_review.pdf"
    page_pngs = sorted((pkg / "qa" / "docx_render").glob("page-*.png"))
    ok = zipfile.is_zipfile(docx) and pdf.read_bytes().startswith(b"%PDF") and len(page_pngs) == 16
    add(checks, "review document files", ok, f"DOCX zip={zipfile.is_zipfile(docx)}; PDF header={pdf.read_bytes().startswith(b'%PDF')}; QA pages={len(page_pngs)}")


def manifest_checksum_checks(pkg: Path, checks: list[Check]) -> None:
    manifest = (pkg / "PACKAGE_MANIFEST.txt").read_text(encoding="utf-8").splitlines()
    actual = sorted(str(path.relative_to(pkg)) for path in pkg.rglob("*") if path.is_file())
    add(checks, "manifest matches files", manifest == actual, f"manifest={len(manifest)} files; actual={len(actual)} files")

    checksum_lines = (pkg / "release" / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    checksum_map = {}
    for line in checksum_lines:
        digest, rel = line.split(maxsplit=1)
        checksum_map[rel] = digest
    failures: list[str] = []
    for rel, expected_digest in checksum_map.items():
        path = pkg / rel
        if not path.exists():
            failures.append(f"missing {rel}")
        elif sha256(path) != expected_digest:
            failures.append(f"checksum mismatch {rel}")
    add(checks, "per-file SHA256", not failures, "failures: " + "; ".join(failures[:10]) if failures else f"{len(checksum_map)} checksums verified")


def write_report(path: Path, zip_path: Path, checks: list[Check]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    failures = [check for check in checks if check.status != "PASS"]
    lines = [
        "# CEE v0.3 Internal Review Package Smoke Test",
        "",
        "Date: 2026-06-05.",
        "",
        f"Package zip: `{display_path(zip_path)}`.",
        f"Package SHA256: `{sha256(zip_path)}`.",
        "",
        "Status: " + ("passed" if not failures else "failed"),
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    for check in checks:
        lines.append(f"| {check.name} | {check.status} | {check.detail} |")
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "This smoke test verifies internal review package integrity and derived-output consistency. It does not rerun the full raw-data MF2013 or Zhao workflows because the raw J-SHIS/NIED flatfile is intentionally excluded from the package.",
            "",
            "## Integrity Note",
            "",
            "Passing this smoke test does not make the manuscript submission-ready. Domain-expert review, reference-manager verification, public repository/DOI setup, and author/legal review remain required.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    zip_path = DEFAULT_ZIP
    extract_root = DEFAULT_EXTRACT_DIR
    report_path = DEFAULT_REPORT

    if extract_root.exists():
        shutil.rmtree(extract_root)
    extract_root.mkdir(parents=True)

    checks: list[Check] = []
    add(checks, "zip exists", zip_path.exists(), display_path(zip_path))
    add(checks, "zip integrity", zipfile.is_zipfile(zip_path), "zipfile header valid")

    with zipfile.ZipFile(zip_path) as archive:
        bad_member = archive.testzip()
        add(checks, "zip member CRC", bad_member is None, "all members OK" if bad_member is None else f"bad member: {bad_member}")
        archive.extractall(extract_root)

    pkg = extract_root / PACKAGE_DIRNAME
    add(checks, "extracted package directory", pkg.exists(), display_path(pkg))

    expected_file_checks(pkg, checks)
    forbidden_scan(pkg, checks)
    script_compile_checks(pkg, checks)
    table_checks(pkg, checks)
    figure_and_document_checks(pkg, checks)
    manifest_checksum_checks(pkg, checks)

    write_report(report_path, zip_path, checks)
    print(report_path)
    failures = [check for check in checks if check.status != "PASS"]
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
