#!/usr/bin/env python3
"""Rule-based AI-term sensitivity for the official J-SHIS/MF2013 residual result.

This is not a full MF2013 reproduction. The public flatfile exposes the
volcanic-front distances needed for the AI term, but not a direct plate-membership
flag. This script therefore applies a transparent geographic rule set only as a
sensitivity check and keeps PH omitted.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ZIP = PROJECT_ROOT / "work" / "external_data" / "jshis_gmf" / "flatfile_sub1-v2024.zip"
DEFAULT_COEFFS = PROJECT_ROOT / "work" / "external_data" / "jshis_mf2013" / "MF13rev_coefs.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"

TARGETS = [
    {
        "target": "maxaccrd050",
        "target_label": "PGA RotD50",
        "period": "Acc.",
        "observed_col": "maxaccrd050",
    },
    {
        "target": "rsaccrd050d005t0030",
        "target_label": "SA(0.3s) RotD50",
        "period": "0.30",
        "observed_col": "rsaccrd050d005t0030",
    },
    {
        "target": "rsaccrd050d005t0100",
        "target_label": "SA(1.0s) RotD50",
        "period": "1.00",
        "observed_col": "rsaccrd050d005t0100",
    },
    {
        "target": "rsaccrd050d005t0300",
        "target_label": "SA(3.0s) RotD50",
        "period": "3.00",
        "observed_col": "rsaccrd050d005t0300",
    },
]

SMREC_USECOLS = [
    "smrec_id",
    "site_id",
    "eq_source_id",
    "maxaccrd050",
    "rsaccrd050d005t0030",
    "rsaccrd050d005t0100",
    "rsaccrd050d005t0300",
    "fault_dist",
]

SITE_USECOLS = [
    "siteid2",
    "lat",
    "lon",
    "avs30",
    "d1400",
    "dist_vf_mf13_nejapan",
    "dist_vf_mf13_swjapan",
]

SOURCE_USECOLS = [
    "eq_source_id",
    "segment_idx",
    "jem_lat",
    "jem_lon",
    "jem_depth",
    "mw",
    "eq_location_type_id",
]

MODELS = [
    ("mf2013_basic", False, False),
    ("mf2013_site", True, False),
    ("mf2013_site_ai_rule_sensitivity", True, True),
]

SOURCE_CLASS_LABELS = {
    1: "CRUSTAL",
    2: "INTERPLATE",
    3: "INTRAPLATE",
}


def load_coeffs(path: Path) -> pd.DataFrame:
    coeffs = pd.read_csv(path)
    coeffs.columns = [col.strip() for col in coeffs.columns]
    for col in coeffs.columns:
        if col != "Period":
            coeffs[col] = pd.to_numeric(coeffs[col], errors="coerce")
    coeffs["Period"] = coeffs["Period"].astype(str).str.strip()
    return coeffs.set_index("Period")


def read_schema_tables(zip_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(zip_path) as zf:
        with zf.open("site_schema.tsv") as handle:
            site = pd.read_csv(handle, sep="\t", usecols=SITE_USECOLS)
        with zf.open("source_schema.tsv") as handle:
            source = pd.read_csv(handle, sep="\t", usecols=SOURCE_USECOLS)
    for col in SITE_USECOLS:
        site[col] = pd.to_numeric(site[col], errors="coerce")
    site["siteid2"] = site["siteid2"].astype("Int64")
    site = site.drop_duplicates("siteid2").copy()
    for col in SOURCE_USECOLS:
        source[col] = pd.to_numeric(source[col], errors="coerce")
    source["eq_source_id"] = source["eq_source_id"].astype("Int64")
    return source.sort_values(["eq_source_id", "segment_idx"]).drop_duplicates("eq_source_id"), site


def source_class_terms(source_class: pd.Series, coeff: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    k = source_class.astype("Int64")
    b = np.select(
        [k.eq(1), k.eq(2), k.eq(3)],
        [coeff["b1"], coeff["b2"], coeff["b3"]],
        default=np.nan,
    )
    c = np.select(
        [k.eq(1), k.eq(2), k.eq(3)],
        [coeff["c1"], coeff["c2"], coeff["c3"]],
        default=np.nan,
    )
    return b, c


def ai_masks(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Transparent approximate applicability masks from the public MF2013 page."""

    depth = pd.to_numeric(frame["jem_depth"], errors="coerce")
    source_lon = pd.to_numeric(frame["jem_lon"], errors="coerce")
    site_lon = pd.to_numeric(frame["lon"], errors="coerce")
    site_lat = pd.to_numeric(frame["lat"], errors="coerce")
    source_class = pd.to_numeric(frame["eq_location_type_id"], errors="coerce")

    subduction_like = source_class.isin([2, 3])
    ne_mask = subduction_like & depth.gt(30.0) & site_lat.ge(36.0) & source_lon.ge(136.9)
    sw_mask = subduction_like & depth.ge(60.0) & site_lon.lt(136.9) & source_lon.lt(136.9)
    return ne_mask.fillna(False), sw_mask.fillna(False)


def ai_term(frame: pd.DataFrame, coeff: pd.Series) -> pd.Series:
    depth = pd.to_numeric(frame["jem_depth"], errors="coerce")
    ne_mask, sw_mask = ai_masks(frame)
    ai = pd.Series(0.0, index=frame.index, dtype=float)
    ne_xvf = pd.to_numeric(frame["dist_vf_mf13_nejapan"], errors="coerce")
    sw_xvf = pd.to_numeric(frame["dist_vf_mf13_swjapan"], errors="coerce").clip(upper=75.0)
    ai.loc[ne_mask] = coeff["gNE"] * ne_xvf.loc[ne_mask] * (depth.loc[ne_mask] - 30.0)
    ai.loc[sw_mask] = coeff["gSW"] * sw_xvf.loc[sw_mask] * (depth.loc[sw_mask] - 30.0)
    return ai


def mf2013_prediction(
    frame: pd.DataFrame,
    coeff: pd.Series,
    include_site: bool,
    include_ai: bool,
) -> pd.Series:
    mw = np.minimum(pd.to_numeric(frame["mw"], errors="coerce").astype(float), 8.2)
    x = np.maximum(pd.to_numeric(frame["fault_dist"], errors="coerce").astype(float), 1.0)
    b, c = source_class_terms(frame["eq_location_type_id"], coeff)
    pred = (
        coeff["a"] * (mw - 16.0) ** 2
        + b * x
        + c
        - np.log10(x + coeff["d"] * 10.0 ** (0.5 * mw))
    )
    pred = pd.Series(pred, index=frame.index)
    if include_site:
        d1400 = pd.to_numeric(frame["d1400"], errors="coerce").astype(float)
        avs30 = pd.to_numeric(frame["avs30"], errors="coerce").astype(float)
        gd = pd.Series(np.nan, index=frame.index, dtype=float)
        d1400_mask = d1400.notna() & d1400.gt(0)
        gd.loc[d1400_mask] = coeff["pd"] * np.log10(
            np.maximum(coeff["Dlmin"], d1400.loc[d1400_mask]) / 300.0
        )
        gs = pd.Series(np.nan, index=frame.index, dtype=float)
        avs30_mask = avs30.notna() & avs30.gt(0)
        gs.loc[avs30_mask] = coeff["ps"] * np.log10(
            np.minimum(coeff["Vsmax"], avs30.loc[avs30_mask]) / 350.0
        )
        pred = pred + gd + gs
    if include_ai:
        pred = pred + ai_term(frame, coeff)
    return pred


def target_model_metrics(
    joined: pd.DataFrame,
    coeffs: pd.DataFrame,
    target_info: dict[str, str],
    model_name: str,
    include_site: bool,
    include_ai: bool,
) -> dict[str, object]:
    coeff = coeffs.loc[target_info["period"]]
    observed_raw = pd.to_numeric(joined[target_info["observed_col"]], errors="coerce")
    observed = np.log10(observed_raw.where(observed_raw > 0))
    predicted = mf2013_prediction(joined, coeff, include_site=include_site, include_ai=include_ai)
    usable = (
        observed.notna()
        & predicted.notna()
        & pd.to_numeric(joined["fault_dist"], errors="coerce").gt(0)
        & pd.to_numeric(joined["mw"], errors="coerce").notna()
        & joined["eq_location_type_id"].isin([1, 2, 3])
    )
    residual = observed.loc[usable] - predicted.loc[usable]
    ne_mask, sw_mask = ai_masks(joined.loc[usable])
    return {
        "target": target_info["target"],
        "target_label": target_info["target_label"],
        "period": target_info["period"],
        "model": model_name,
        "include_site": include_site,
        "include_ai_rule_sensitivity": include_ai,
        "n_records": int(usable.sum()),
        "residual_sum": float(residual.sum()),
        "residual_sq_sum": float((residual**2).sum()),
        "absolute_error_sum": float(residual.abs().sum()),
        "ai_ne_rule_records": int(ne_mask.sum()) if include_ai else 0,
        "ai_sw_rule_records": int(sw_mask.sum()) if include_ai else 0,
    }


def finalize_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    n = out["n_records"].replace(0, np.nan)
    out["mean_residual"] = out["residual_sum"] / n
    out["mae"] = out["absolute_error_sum"] / n
    out["rmse"] = np.sqrt(out["residual_sq_sum"] / n)
    rows = []
    for target, sub in out.groupby("target", dropna=False):
        basic = sub[sub["model"].eq("mf2013_basic")]
        site = sub[sub["model"].eq("mf2013_site")]
        if basic.empty:
            continue
        basic_mae = float(basic.iloc[0]["mae"])
        site_mae = float(site.iloc[0]["mae"]) if not site.empty else np.nan
        for _, row in sub.iterrows():
            item = row.to_dict()
            item["mae_reduction_vs_basic"] = basic_mae - row["mae"]
            item["percent_mae_reduction_vs_basic"] = (basic_mae - row["mae"]) / basic_mae * 100.0
            item["mae_reduction_vs_site"] = site_mae - row["mae"] if np.isfinite(site_mae) else np.nan
            item["percent_mae_reduction_vs_site"] = (
                (site_mae - row["mae"]) / site_mae * 100.0 if np.isfinite(site_mae) else np.nan
            )
            rows.append(item)
    return pd.DataFrame(rows)


def write_report(metrics: pd.DataFrame, output_dir: Path) -> None:
    lines = ["# MF2013 AI Rule-Sensitivity Audit\n\n"]
    lines.append("## Scope\n")
    lines.append("- This is a sensitivity check, not a full MF2013 or full J-SHIS hazard-map reproduction.\n")
    lines.append("- Models compared: `mf2013_basic`, `mf2013_site`, and `mf2013_site_ai_rule_sensitivity`.\n")
    lines.append("- PH is not applied because the public schema used here does not expose a direct Philippine Sea Plate membership flag.\n\n")
    lines.append("## AI rule used\n")
    lines.append("- NE rule: subduction-like source class, hypocenter depth > 30 km, site latitude >= 36.0, source longitude >= 136.9.\n")
    lines.append("- SW rule: subduction-like source class, hypocenter depth >= 60 km, site longitude < 136.9, source longitude < 136.9.\n")
    lines.append("- SW volcanic-front distance is capped at 75 km, following the public MF2013 page statement for southwestern Japan.\n")
    lines.append("- These rules approximate public applicability language and should remain labelled as sensitivity rules.\n\n")

    lines.append("## Overall metrics\n")
    for _, row in metrics.sort_values(["target_label", "model"]).iterrows():
        lines.append(
            f"- {row['target_label']} / {row['model']}: n={int(row['n_records']):,}, "
            f"MAE={row['mae']:.4f}, bias={row['mean_residual']:.4f}, "
            f"MAE reduction vs basic={row['percent_mae_reduction_vs_basic']:.1f}%, "
            f"MAE reduction vs site={row['percent_mae_reduction_vs_site']:.2f}%.\n"
        )

    ai = metrics[metrics["model"].eq("mf2013_site_ai_rule_sensitivity")]
    lines.append("\n## Affected records\n")
    for _, row in ai.sort_values("target_label").iterrows():
        lines.append(
            f"- {row['target_label']}: NE-rule records={int(row['ai_ne_rule_records']):,}, "
            f"SW-rule records={int(row['ai_sw_rule_records']):,}.\n"
        )

    lines.append("\n## Interpretation\n")
    site = metrics[metrics["model"].eq("mf2013_site")]
    site_ai = metrics[metrics["model"].eq("mf2013_site_ai_rule_sensitivity")]
    merged = site.merge(site_ai, on=["target", "target_label"], suffixes=("_site", "_ai"))
    for _, row in merged.sort_values("target_label").iterrows():
        delta = row["percent_mae_reduction_vs_basic_ai"] - row["percent_mae_reduction_vs_basic_site"]
        lines.append(
            f"- {row['target_label']}: adding the AI sensitivity changes percent MAE reduction by {delta:.2f} percentage points.\n"
        )
    lines.append(
        "\nThe acceptance criterion for the CEE storyline is whether SA(1.0 s) and SA(3.0 s) remain much stronger than PGA after this sensitivity. "
        "This audit should not be used to claim official PH-corrected MF2013 reproduction.\n"
    )
    lines.append("\n## Outputs\n")
    lines.append("- `jshis_mf2013_ai_sensitivity_metrics.csv`\n")
    lines.append("- `jshis_mf2013_ai_sensitivity.md`\n")
    (output_dir / "jshis_mf2013_ai_sensitivity.md").write_text("".join(lines), encoding="utf-8")


def run(zip_path: Path, coeff_path: Path, output_dir: Path, chunksize: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    coeffs = load_coeffs(coeff_path)
    source, site = read_schema_tables(zip_path)
    pieces: list[pd.DataFrame] = []

    with ZipFile(zip_path) as zf:
        with zf.open("smrec_schema.tsv") as handle:
            reader = pd.read_csv(
                handle,
                sep="\t",
                usecols=SMREC_USECOLS,
                chunksize=chunksize,
                low_memory=False,
            )
            for chunk_index, chunk in enumerate(reader, start=1):
                chunk["site_id"] = pd.to_numeric(chunk["site_id"], errors="coerce").astype("Int64")
                chunk["siteid2"] = (chunk["site_id"] // 10).astype("Int64")
                chunk["eq_source_id"] = pd.to_numeric(chunk["eq_source_id"], errors="coerce").astype("Int64")
                joined = chunk.merge(site, on="siteid2", how="left").merge(source, on="eq_source_id", how="left")
                rows = []
                for target_info in TARGETS:
                    for model_name, include_site, include_ai in MODELS:
                        rows.append(
                            target_model_metrics(
                                joined,
                                coeffs,
                                target_info,
                                model_name,
                                include_site,
                                include_ai,
                            )
                        )
                pieces.append(pd.DataFrame(rows))
                print(f"chunk={chunk_index} joined_rows={len(joined):,}", flush=True)

    raw = pd.concat(pieces, ignore_index=True)
    summed = (
        raw.groupby(
            ["target", "target_label", "period", "model", "include_site", "include_ai_rule_sensitivity"],
            as_index=False,
        )[
            [
                "n_records",
                "residual_sum",
                "residual_sq_sum",
                "absolute_error_sum",
                "ai_ne_rule_records",
                "ai_sw_rule_records",
            ]
        ]
        .sum()
    )
    metrics = finalize_metrics(summed)
    metrics.to_csv(output_dir / "jshis_mf2013_ai_sensitivity_metrics.csv", index=False)
    write_report(metrics, output_dir)
    print(f"metric_rows={len(metrics):,}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--coeffs", type=Path, default=DEFAULT_COEFFS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunksize", type=int, default=100_000)
    args = parser.parse_args()
    run(args.zip, args.coeffs, args.output_dir, args.chunksize)


if __name__ == "__main__":
    main()
