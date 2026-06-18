#!/usr/bin/env python3
"""Audit consistency between K-NET header PGA and waveform-derived peaks."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_record_waveform_features.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def rel_diff(a: pd.Series, b: pd.Series) -> pd.Series:
    return (a - b).abs() / b.replace(0, np.nan).abs()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures = args.output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, dtype={"_id": str, "EventName": str})
    comp_cols = [
        "pga_component_gal_NS",
        "pga_component_gal_EW",
        "pga_component_gal_UD",
        "pga_horizontal_from_waveform_gal",
        "pga_3c_from_waveform_gal",
    ]
    labels = ["NS", "EW", "UD", "H", "3C"]
    diffs = pd.concat([rel_diff(df[col], df["PGA_gal"]) for col in comp_cols], axis=1)
    diffs.columns = labels
    df["closest_pga_definition"] = diffs.idxmin(axis=1)
    df["closest_pga_relative_diff"] = diffs.min(axis=1)
    df["header_matches_horizontal_1pct"] = rel_diff(
        df["pga_horizontal_from_waveform_gal"], df["PGA_gal"]
    ) <= 0.01
    df["header_matches_3c_1pct"] = rel_diff(
        df["pga_3c_from_waveform_gal"], df["PGA_gal"]
    ) <= 0.01
    df["header_matches_any_component_1pct"] = df["closest_pga_relative_diff"] <= 0.01
    df["horizontal_header_ratio"] = df["pga_horizontal_from_waveform_gal"] / df["PGA_gal"]
    df["three_component_header_ratio"] = df["pga_3c_from_waveform_gal"] / df["PGA_gal"]

    audit_cols = [
        "_id",
        "EventName",
        "StationCode",
        "Magnitude",
        "Distance",
        "PGA_gal",
        "pga_component_gal_NS",
        "pga_component_gal_EW",
        "pga_component_gal_UD",
        "pga_horizontal_from_waveform_gal",
        "pga_3c_from_waveform_gal",
        "horizontal_header_ratio",
        "three_component_header_ratio",
        "closest_pga_definition",
        "closest_pga_relative_diff",
        "header_matches_horizontal_1pct",
        "header_matches_3c_1pct",
        "header_matches_any_component_1pct",
    ]
    df[audit_cols].to_csv(args.output_dir / "knet_pga_consistency_audit.csv", index=False)

    plt.figure(figsize=(7.2, 4.5))
    ratio = df["horizontal_header_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
    plt.hist(np.clip(ratio, 0, 3), bins=80, color="#4c78a8")
    plt.axvline(1.0, color="k", linewidth=1)
    plt.xlabel("Waveform horizontal PGA / header PGA_gal")
    plt.ylabel("Record count")
    plt.title("K-NET header PGA consistency")
    plt.tight_layout()
    plt.savefig(figures / "knet_header_waveform_pga_ratio_histogram.png", dpi=220)
    plt.close()

    counts = df["closest_pga_definition"].value_counts().reindex(labels, fill_value=0)
    plt.figure(figsize=(6.8, 4.2))
    plt.bar(counts.index, counts.values, color="#59a14f")
    plt.ylabel("Record count")
    plt.xlabel("Closest waveform definition to header PGA")
    plt.title("Which waveform peak matches header PGA?")
    plt.tight_layout()
    plt.savefig(figures / "knet_header_pga_closest_definition.png", dpi=220)
    plt.close()

    lines = ["# K-NET Header PGA Consistency Audit\n\n"]
    lines.append("## Summary\n")
    lines.append(f"- Records audited: {len(df):,}\n")
    lines.append(
        f"- Header matches waveform horizontal PGA within 1%: {int(df['header_matches_horizontal_1pct'].sum()):,} / {len(df):,}\n"
    )
    lines.append(
        f"- Header matches waveform 3C PGA within 1%: {int(df['header_matches_3c_1pct'].sum()):,} / {len(df):,}\n"
    )
    lines.append(
        f"- Header matches any single/combined waveform PGA definition within 1%: {int(df['header_matches_any_component_1pct'].sum()):,} / {len(df):,}\n"
    )
    lines.append("- Closest definition counts:\n")
    for key, value in counts.items():
        lines.append(f"  - {key}: {int(value):,}\n")
    lines.append("\n## Ratio distribution\n")
    lines.append(f"- Median horizontal/header ratio: {ratio.median():.6f}\n")
    lines.append(f"- 2.5%-97.5% range: {ratio.quantile(0.025):.6f} to {ratio.quantile(0.975):.6f}\n")
    lines.append(f"- 1%-99% range: {ratio.quantile(0.01):.6f} to {ratio.quantile(0.99):.6f}\n")
    lines.append("\n## Recommendation\n")
    lines.append(
        "- Use waveform-derived horizontal, vertical, and 3C peaks as primary targets/features in future models.\n"
    )
    lines.append(
        "- Treat header `PGA_gal` as a useful metadata field but not a uniform definition across all records.\n"
    )
    lines.append(
        "- Keep records with large discrepancies flagged for quality control instead of silently discarding them.\n"
    )
    (args.output_dir / "knet_pga_consistency_audit.md").write_text(
        "".join(lines), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
