#!/usr/bin/env python3
"""Audit J-SHIS flatfile integrity and the MF2013 record-selection flow."""

from __future__ import annotations

import argparse
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as model


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT = ARTICLE / "supplement"
DEFAULT_FLATFILE = ROOT / "work" / "external_data" / "jshis_gmf" / "flatfile_sub1-v2024.zip"
DEFAULT_COEFFICIENTS = ROOT / "work" / "external_data" / "jshis_mf2013" / "MF13rev_coefs.csv"
DEFAULT_OUTPUT_CSV = SUPPLEMENT / "jshis_flatfile_selection_audit.csv"
DEFAULT_OUTPUT_REPORT = SUPPLEMENT / "jshis_flatfile_selection_audit.md"


def mask_summary(frame: pd.DataFrame, mask: pd.Series) -> tuple[int, int, int]:
    selected = frame.loc[mask]
    return len(selected), selected["eq_source_id"].nunique(), selected["siteid2"].nunique()


def add_stage(
    rows: list[dict[str, object]],
    frame: pd.DataFrame,
    stage_order: int,
    stage: str,
    mask: pd.Series,
    previous_count: int | None,
) -> int:
    records, earthquakes, stations = mask_summary(frame, mask)
    rows.append(
        {
            "stage_order": stage_order,
            "stage": stage,
            "records": records,
            "excluded_at_stage": 0 if previous_count is None else previous_count - records,
            "earthquakes": earthquakes,
            "stations": stations,
        }
    )
    return records


def run(args: argparse.Namespace) -> None:
    coefficients = model.load_coefficients(args.coefficients)
    site, source = model.load_metadata(args.flatfile)
    response_columns = [spec.rotd100_col for spec in model.PERIODS]

    with zipfile.ZipFile(args.flatfile) as archive:
        source_jma = pd.read_csv(
            archive.open("source_schema.tsv"),
            sep="\t",
            usecols=["eq_source_id", "segment_idx", "mjma"],
        )
        records = pd.read_csv(
            archive.open("smrec_schema.tsv"),
            sep="\t",
            usecols=["smrec_id", "site_id", "eq_source_id", "fault_dist", *response_columns],
            low_memory=False,
        )

    for column in ["eq_source_id", "segment_idx", "mjma"]:
        source_jma[column] = pd.to_numeric(source_jma[column], errors="coerce")
    source_jma = (
        source_jma.sort_values(["eq_source_id", "segment_idx"])
        .drop_duplicates("eq_source_id", keep="first")
        .drop(columns="segment_idx")
    )
    source = source.merge(source_jma, on="eq_source_id", how="left", validate="one_to_one")

    records["site_id"] = pd.to_numeric(records["site_id"], errors="coerce").astype("Int64")
    records["siteid2"] = (records["site_id"] // 10).astype("Int64")
    records["eq_source_id"] = pd.to_numeric(records["eq_source_id"], errors="coerce").astype("Int64")

    if records["smrec_id"].isna().any() or records["smrec_id"].duplicated().any():
        raise ValueError("smrec_id must be complete and unique")

    joined = records.merge(site, on="siteid2", how="left", validate="many_to_one", indicator="_site_join")
    missing_site = int(joined["_site_join"].ne("both").sum())
    joined = joined.drop(columns="_site_join").merge(
        source,
        on="eq_source_id",
        how="left",
        validate="many_to_one",
        indicator="_source_join",
    )
    missing_source = int(joined["_source_join"].ne("both").sum())
    joined = joined.drop(columns="_source_join")
    if missing_site or missing_source:
        raise ValueError(f"incomplete joins: site={missing_site}, source={missing_source}")

    all_records = pd.Series(True, index=joined.index)
    surface = joined["installation_situation_id"].eq(model.PRIMARY_INSTALLATION_SITUATION_ID)
    finite_mw = np.isfinite(pd.to_numeric(joined["mw"], errors="coerce"))
    source_and_distance = joined["eq_location_type_id"].isin([1, 2, 3]) & pd.to_numeric(
        joined["fault_dist"], errors="coerce"
    ).gt(0)
    positive_avs30 = pd.to_numeric(joined["avs30"], errors="coerce").gt(0)

    stage_rows: list[dict[str, object]] = []
    previous = add_stage(stage_rows, joined, 1, "Public sub1-v2024 archive", all_records, None)
    mask = surface.copy()
    previous = add_stage(stage_rows, joined, 2, "Ground-surface installation", mask, previous)
    mask &= finite_mw
    previous = add_stage(stage_rows, joined, 3, "Finite F-net moment magnitude Mw", mask, previous)
    mask &= source_and_distance
    previous = add_stage(stage_rows, joined, 4, "Supported source class and positive distance", mask, previous)
    mask &= positive_avs30
    previous = add_stage(stage_rows, joined, 5, "Positive AVS30", mask, previous)

    period_masks: dict[float, pd.Series] = {}
    for spec in model.PERIODS:
        observed = pd.to_numeric(joined[spec.rotd100_col], errors="coerce")
        basic, site_prediction, site_ai_prediction = model.mf2013_predictions(
            joined, coefficients.loc[spec.coeff_key]
        )
        period_masks[spec.period_s] = (
            mask
            & observed.gt(0)
            & np.isfinite(observed)
            & np.isfinite(basic)
            & np.isfinite(site_prediction)
            & np.isfinite(site_ai_prediction)
        )
    complete_periods = pd.concat(period_masks, axis=1).all(axis=1)
    add_stage(
        stage_rows,
        joined,
        6,
        "Finite RotD100 and MF2013 predictions at all eight periods",
        complete_periods,
        previous,
    )

    period_counts = {period: int(period_mask.sum()) for period, period_mask in period_masks.items()}
    if len(set(period_counts.values())) != 1:
        raise ValueError(f"period-specific sample sizes differ: {period_counts}")

    stage_frame = pd.DataFrame(stage_rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    stage_frame.to_csv(args.output_csv, index=False)

    mjma = pd.to_numeric(joined["mjma"], errors="coerce")
    fault_distance = pd.to_numeric(joined["fault_dist"], errors="coerce")
    surface_nonpositive_avs30 = surface & ~positive_avs30
    overlap = surface_nonpositive_avs30 & ~finite_mw
    lines = [
        "# J-SHIS flatfile selection audit",
        "",
        f"- Archive: `{args.flatfile.relative_to(ROOT)}`",
        f"- Record rows and unique `smrec_id` values: {len(joined):,}",
        f"- Missing site joins: {missing_site:,}",
        f"- Missing source joins: {missing_source:,}",
        f"- JMA magnitude range: {mjma.min():.1f}--{mjma.max():.1f}",
        f"- Shortest-fault-distance range: {fault_distance.min():.4f}--{fault_distance.max():.4f} km",
        f"- Surface records with nonpositive AVS30: {int(surface_nonpositive_avs30.sum()):,}",
        f"- Nonpositive-AVS30 records that also lack finite Mw: {int(overlap.sum()):,}",
        "",
        "## Sequential selection",
        "",
        "| Stage | Records | Excluded at stage | Earthquakes | Stations |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in stage_rows:
        lines.append(
            f"| {row['stage']} | {row['records']:,} | {row['excluded_at_stage']:,} | "
            f"{row['earthquakes']:,} | {row['stations']:,} |"
        )
    lines.extend(
        [
            "",
            "## Period-specific final counts",
            "",
            ", ".join(f"{period:g} s: {count:,}" for period, count in period_counts.items()),
            "",
            "Every record in the downloaded subset satisfies `mjma >= 5`. MF2013 uses the separate F-net `mw` "
            "field, which creates the magnitude-related loss between the surface archive and the model sample.",
        ]
    )
    args.output_report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(stage_frame.to_string(index=False))
    print(f"wrote {args.output_csv.relative_to(ROOT)}")
    print(f"wrote {args.output_report.relative_to(ROOT)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=DEFAULT_COEFFICIENTS)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT_CSV)
    parser.add_argument("--output-report", type=Path, default=DEFAULT_OUTPUT_REPORT)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
