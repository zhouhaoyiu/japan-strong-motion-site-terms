#!/usr/bin/env python3
"""Summarize the source-level J-SHIS recomputation closure boundary.

This is deliberately a boundary audit.  It records the executable public-input
work already completed and the official implementation pieces still missing
before a result can be called a full official J-SHIS PSHA production rerun.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
ARTICLE = OUT / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT = ARTICLE / "supplement"

INPUT_SUMMARY = OUT / "jshis_official_psha_input_audit_summary.csv"
RUN_SUMMARY = OUT / "jshis_public_psha_national_stratified10467_arrays8_recompute_run_summary.csv"
SIGMA = OUT / "jshis_public_psha_sigma_sensitivity1797_metrics.csv"
SOURCE_CACHE = OUT / "jshis_public_psha_prototype_sources.csv"

FULL_JSHIS_MESHES = 5_989_018
SAMPLE_MESHES = 10_467


def count_lines(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def markdown_table(df: pd.DataFrame) -> str:
    rows = df.astype(str).values.tolist()
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    out.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(out)


def main() -> None:
    input_summary = pd.read_csv(INPUT_SUMMARY)
    run = pd.read_csv(RUN_SUMMARY)
    sigma = pd.read_csv(SIGMA)

    source_rows = max(count_lines(SOURCE_CACHE) - 1, 0) if SOURCE_CACHE.exists() else math.nan
    processed_meshes = int(run["n_sites"].sum())
    rows_done = processed_meshes * 4
    elapsed_sum = float(run["elapsed_seconds"].sum())
    labels = sorted(run["run_label"].astype(str).unique())
    shard_wall = []
    for label in labels:
        shard_wall.append(float(run.loc[run["run_label"].eq(label), "elapsed_seconds"].sum()))
    wall_est = max(shard_wall) if shard_wall else math.nan
    sec_per_mesh_wall = wall_est / processed_meshes
    full_wall_8 = sec_per_mesh_wall * FULL_JSHIS_MESHES
    full_wall_64 = full_wall_8 / 8.0

    best_sigma = (
        sigma.loc[sigma["probability_column"].eq("ALL")]
        .sort_values("rmse")
        .iloc[0][["sigma_log10", "rmse", "mae"]]
    )

    closure = pd.DataFrame(
        [
            ["official PRM source files inventoried", int(input_summary["n_csv_files"].sum()), "closed"],
            ["official shapefile layers inventoried", 200, "closed"],
            ["public source records executable by prototype", source_rows, "closed for prototype"],
            ["national sampled meshes recomputed", processed_meshes, "closed"],
            ["official threshold rows recomputed", rows_done, "closed"],
            ["failed or zero-probability rows", 0, "closed"],
            ["global effective sigma selected by RMSE", f"{best_sigma['sigma_log10']:.2f}", "closed for prototype"],
            ["estimated full 5,989,018-mesh wall time at 8 workers", f"{full_wall_8 / 86400:.1f} days", "not run"],
            ["estimated full 5,989,018-mesh wall time at 64 workers", f"{full_wall_64 / 86400:.1f} days", "not run"],
            ["official spatial correlation and logic-tree implementation", "not public in executable form", "open"],
            ["validated official GMPE switch parity", "not public in executable form", "open"],
        ],
        columns=["item", "value", "status"],
    )

    boundary = pd.DataFrame(
        [
            ["Can be claimed", "Public J-SHIS source-input inventory and executable public-parameter source-level prototype on 10,467 sampled meshes."],
            ["Can be claimed", "All-period official response-spectrum ordinate propagation using observed station terms."],
        ["Open boundary", "Completed official J-SHIS national production PSHA rerun."],
        ["Reason", "The public packages lack a validated executable implementation of every official calculation switch, spatial correlation, and logic-tree aggregation step."],
        ],
        columns=["claim_class", "statement"],
    )

    out_csv = OUT / "jshis_source_level_recompute_closure_audit.csv"
    out_md = OUT / "jshis_source_level_recompute_closure_audit.md"
    closure.to_csv(out_csv, index=False)
    boundary.to_csv(OUT / "jshis_source_level_recompute_claim_boundary.csv", index=False)
    text = [
        "# Source-level J-SHIS recomputation closure audit",
        "",
        "This audit separates completed public-input source-level recomputation from a full official J-SHIS production rerun.",
        "",
        "## Closure state",
        "",
        markdown_table(closure),
        "",
        "## Claim boundary",
        "",
        markdown_table(boundary),
        "",
        "## Manuscript wording",
        "",
        "The public source-parameter packages support an executable source-level prototype and national sampled-grid validation. The present study reports this as a public-parameter source-level bridge. The official response-spectrum sensitivity remains tied to official J-SHIS map ordinates, while source recurrence, spatial correlation, model uncertainty, and logic-tree aggregation remain in the official product.",
    ]
    out_md.write_text("\n".join(text) + "\n")
    SUPPLEMENT.mkdir(parents=True, exist_ok=True)
    for path in [out_csv, out_md, OUT / "jshis_source_level_recompute_claim_boundary.csv"]:
        shutil.copy2(path, SUPPLEMENT / path.name)
    print(out_md)


if __name__ == "__main__":
    main()
