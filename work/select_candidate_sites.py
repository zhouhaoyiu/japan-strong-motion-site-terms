#!/usr/bin/env python3
"""Select candidate stations for the source-path-site/site-response study."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import pandas as pd


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def code_family(code: str) -> str:
    code = str(code)
    return code[:3]


def station_number(code: str) -> str:
    digits = "".join(ch for ch in str(code) if ch.isdigit())
    return digits[-2:] if len(digits) >= 2 else digits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    knet = pd.read_csv(args.output_dir / "knet_station_residual_terms.csv")
    kik = pd.read_csv(args.output_dir / "kiknet_surface_downhole_ratios.csv")

    knet["code_family"] = knet["StationCode"].map(code_family)
    knet["station_number_hint"] = knet["StationCode"].map(station_number)
    kik["code_family"] = kik["station"].map(code_family)
    kik["station_number_hint"] = kik["station"].map(station_number)

    knet_supported = knet[knet["n_records"] >= 30].copy()
    knet_high = knet_supported.sort_values("mean_log10_pga_residual", ascending=False).head(40)
    knet_low = knet_supported.sort_values("mean_log10_pga_residual", ascending=True).head(40)

    kik_high = kik.sort_values("horizontal_pga_ratio", ascending=False).head(40)
    kik_low = kik.sort_values("horizontal_pga_ratio", ascending=True).head(40)

    family_summary = (
        knet_supported.groupby("code_family")
        .agg(
            knet_station_count=("StationCode", "nunique"),
            knet_record_count=("n_records", "sum"),
            mean_knet_residual=("mean_log10_pga_residual", "mean"),
            max_knet_residual=("mean_log10_pga_residual", "max"),
        )
        .join(
            kik.groupby("code_family").agg(
                kik_station_count=("station", "nunique"),
                median_kik_pga_ratio=("horizontal_pga_ratio", "median"),
                max_kik_pga_ratio=("horizontal_pga_ratio", "max"),
                median_kik_low_freq_ratio=("spectral_ratio_0.2_0.5_hz", "median"),
            ),
            how="outer",
        )
        .reset_index()
    )
    family_summary["has_knet_and_kiknet"] = (
        family_summary["knet_station_count"].fillna(0).gt(0)
        & family_summary["kik_station_count"].fillna(0).gt(0)
    )
    family_summary = family_summary.sort_values(
        ["has_knet_and_kiknet", "max_knet_residual", "max_kik_pga_ratio"],
        ascending=[False, False, False],
    )

    candidates = {
        "knet_high_positive_residual_supported": knet_high,
        "knet_high_negative_residual_supported": knet_low,
        "kiknet_high_surface_downhole_ratio": kik_high,
        "kiknet_low_surface_downhole_ratio": kik_low,
        "regional_family_summary": family_summary,
    }
    for name, frame in candidates.items():
        frame.to_csv(args.output_dir / f"{name}.csv", index=False)

    if importlib.util.find_spec("openpyxl") is not None or importlib.util.find_spec("xlsxwriter") is not None:
        with pd.ExcelWriter(args.output_dir / "candidate_site_tables.xlsx") as writer:
            for sheet, frame in candidates.items():
                frame.to_excel(writer, sheet_name=sheet[:31], index=False)

    knet_high.to_csv(args.output_dir / "candidate_knet_high_residual_stations.csv", index=False)
    kik_high.to_csv(args.output_dir / "candidate_kiknet_high_ratio_stations.csv", index=False)
    family_summary.to_csv(args.output_dir / "candidate_regional_family_summary.csv", index=False)

    lines = ["# Candidate Site Selection\n\n"]
    lines.append("## Purpose\n")
    lines.append(
        "This file prioritizes stations/regions for the next full KiK-net multi-event acquisition. "
        "K-NET station residuals and KiK-net surface/downhole ratios are not treated as the same station unless official metadata later confirm co-location.\n\n"
    )
    lines.append("## Top K-NET stations with positive residuals, n_records >= 30\n")
    for _, row in knet_high.head(12).iterrows():
        lines.append(
            f"- {row['StationCode']}: mean residual={row['mean_log10_pga_residual']:.3f}, "
            f"records={int(row['n_records'])}, family={row['code_family']}\n"
        )
    lines.append("\n## Top KiK-net stations by surface/downhole horizontal PGA ratio\n")
    for _, row in kik_high.head(12).iterrows():
        lines.append(
            f"- {row['station']}: PGA ratio={row['horizontal_pga_ratio']:.2f}, "
            f"RMS ratio={row['horizontal_rms_ratio']:.2f}, family={row['code_family']}\n"
        )
    both = family_summary[family_summary["has_knet_and_kiknet"]].head(12)
    lines.append("\n## Priority code families appearing in both K-NET and local KiK-net data\n")
    for _, row in both.iterrows():
        lines.append(
            f"- {row['code_family']}: K-NET stations={int(row['knet_station_count'])}, "
            f"KiK-net stations={int(row['kik_station_count'])}, "
            f"max K-NET residual={row['max_knet_residual']:.3f}, "
            f"max KiK ratio={row['max_kik_pga_ratio']:.2f}\n"
        )
    lines.append("\n## Guardrail\n")
    lines.append(
        "Code-family matching is only a prioritization heuristic. Publication claims need official station coordinates, borehole depth, and site profiles.\n"
    )
    (args.output_dir / "candidate_site_selection.md").write_text("".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
