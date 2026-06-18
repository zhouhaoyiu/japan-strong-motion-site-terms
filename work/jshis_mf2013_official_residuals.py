#!/usr/bin/env python3
"""Official J-SHIS/NIED flatfile residuals against MF2013 core equations.

This script uses the J-SHIS ground-motion flatfile schemas directly:
`smrec_schema.tsv`, `source_schema.tsv`, and `site_schema.tsv`. It computes raw
residuals against the Morikawa-Fujiwara 2013 coefficient CSV for the core
magnitude-distance equation and for the same equation plus D1400/AVS30 site
terms. It intentionally does not apply AI or PH terms because those require
additional plate/region rules that should be audited separately.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats


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
    "multiple",
]

SITE_USECOLS = [
    "siteid2",
    "site_code",
    "site_name",
    "obs_network_id",
    "vs10",
    "vs20",
    "vs30",
    "avs30",
    "d1100",
    "d1400",
    "d1700",
    "d2100",
    "dbase",
]

SOURCE_USECOLS = [
    "eq_source_id",
    "segment_idx",
    "jem_origin_time",
    "jem_depth",
    "mjma",
    "mw",
    "eq_location_type_id",
    "eq_location_type_id_source",
    "eq_event_name",
]

SOURCE_CLASS_LABELS = {
    1: "CRUSTAL",
    2: "INTERPLATE",
    3: "INTRAPLATE",
    10: "OTHER",
    99: "UNKNOWN",
}

NETWORK_LABELS = {
    1: "K-NET",
    2: "KiK-net",
}

AGG_VALUE_COLS = [
    "n_records",
    "residual_sum",
    "residual_sq_sum",
    "absolute_error_sum",
    "observed_sum",
    "predicted_sum",
]


def load_coeffs(path: Path) -> pd.DataFrame:
    coeffs = pd.read_csv(path)
    coeffs.columns = [col.strip() for col in coeffs.columns]
    coeffs["Period"] = coeffs["Period"].astype(str).str.strip()
    for col in coeffs.columns:
        if col != "Period":
            coeffs[col] = pd.to_numeric(coeffs[col], errors="coerce")
    return coeffs.set_index("Period")


def read_schema_tables(zip_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(zip_path) as zf:
        with zf.open("site_schema.tsv") as handle:
            site = pd.read_csv(handle, sep="\t", usecols=SITE_USECOLS)
        with zf.open("source_schema.tsv") as handle:
            source = pd.read_csv(handle, sep="\t", usecols=SOURCE_USECOLS)

    for col in SITE_USECOLS:
        if col not in ["site_code", "site_name"]:
            site[col] = pd.to_numeric(site[col], errors="coerce")
    site["siteid2"] = site["siteid2"].astype("Int64")
    site = site.drop_duplicates("siteid2").copy()

    for col in SOURCE_USECOLS:
        if col not in ["jem_origin_time", "eq_event_name"]:
            source[col] = pd.to_numeric(source[col], errors="coerce")
    source["eq_source_id"] = source["eq_source_id"].astype("Int64")
    source = source.sort_values(["eq_source_id", "segment_idx"]).drop_duplicates("eq_source_id")
    return site, source


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


def mf2013_prediction(frame: pd.DataFrame, coeff: pd.Series, include_site: bool) -> pd.Series:
    mw = np.minimum(pd.to_numeric(frame["mw"], errors="coerce").astype(float), 8.2)
    x = np.maximum(pd.to_numeric(frame["fault_dist"], errors="coerce").astype(float), 1.0)
    b, c = source_class_terms(frame["eq_location_type_id"], coeff)
    pred = (
        coeff["a"] * (mw - 16.0) ** 2
        + b * x
        + c
        - np.log10(x + coeff["d"] * 10.0 ** (0.5 * mw))
    )
    if include_site:
        d1400 = pd.to_numeric(frame["d1400"], errors="coerce").astype(float)
        avs30 = pd.to_numeric(frame["avs30"], errors="coerce").astype(float)
        gd = pd.Series(np.nan, index=frame.index, dtype=float)
        d1400_mask = d1400.notna() & (d1400 > 0)
        gd.loc[d1400_mask] = coeff["pd"] * np.log10(
            np.maximum(coeff["Dlmin"], d1400.loc[d1400_mask]) / 300.0
        )
        gs = pd.Series(np.nan, index=frame.index, dtype=float)
        avs30_mask = avs30.notna() & (avs30 > 0)
        gs.loc[avs30_mask] = coeff["ps"] * np.log10(
            np.minimum(coeff["Vsmax"], avs30.loc[avs30_mask]) / 350.0
        )
        pred = pred + gd + gs
    return pd.Series(pred, index=frame.index)


def residual_frame(
    joined: pd.DataFrame,
    coeffs: pd.DataFrame,
    target_info: dict[str, str],
    model_name: str,
    include_site: bool,
) -> pd.DataFrame:
    coeff = coeffs.loc[target_info["period"]]
    observed_raw = pd.to_numeric(joined[target_info["observed_col"]], errors="coerce")
    observed = np.log10(observed_raw.where(observed_raw > 0))
    predicted = mf2013_prediction(joined, coeff, include_site=include_site)
    usable = (
        observed.notna()
        & predicted.notna()
        & pd.to_numeric(joined["fault_dist"], errors="coerce").gt(0)
        & pd.to_numeric(joined["mw"], errors="coerce").notna()
        & joined["eq_location_type_id"].isin([1, 2, 3])
    )
    if not usable.any():
        return pd.DataFrame()
    out = joined.loc[
        usable,
        [
            "smrec_id",
            "siteid2",
            "site_code",
            "obs_network_id",
            "eq_source_id",
            "eq_location_type_id",
            "mw",
            "mjma",
            "fault_dist",
        ],
    ].copy()
    out["target"] = target_info["target"]
    out["target_label"] = target_info["target_label"]
    out["period"] = target_info["period"]
    out["model"] = model_name
    out["observed_log10"] = observed.loc[usable].to_numpy()
    out["predicted_log10"] = predicted.loc[usable].to_numpy()
    out["residual_log10"] = out["observed_log10"] - out["predicted_log10"]
    out["absolute_error"] = out["residual_log10"].abs()
    out["squared_error"] = out["residual_log10"] ** 2
    return out


def aggregate_residuals(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=group_cols + AGG_VALUE_COLS)
    grouped = (
        frame.groupby(group_cols, dropna=False)
        .agg(
            n_records=("residual_log10", "size"),
            residual_sum=("residual_log10", "sum"),
            residual_sq_sum=("squared_error", "sum"),
            absolute_error_sum=("absolute_error", "sum"),
            observed_sum=("observed_log10", "sum"),
            predicted_sum=("predicted_log10", "sum"),
        )
        .reset_index()
    )
    return grouped


def combine_aggs(parts: list[pd.DataFrame], group_cols: list[str]) -> pd.DataFrame:
    if not parts:
        return pd.DataFrame(columns=group_cols + AGG_VALUE_COLS)
    combined = pd.concat(parts, ignore_index=True)
    combined = (
        combined.groupby(group_cols, dropna=False)[AGG_VALUE_COLS]
        .sum()
        .reset_index()
    )
    return finalize_metrics(combined)


def finalize_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    n = out["n_records"].replace(0, np.nan)
    out["mean_residual"] = out["residual_sum"] / n
    out["mae"] = out["absolute_error_sum"] / n
    out["rmse"] = np.sqrt(out["residual_sq_sum"] / n)
    out["mean_observed_log10"] = out["observed_sum"] / n
    out["mean_predicted_log10"] = out["predicted_sum"] / n
    denom = (out["n_records"] - 1).replace(0, np.nan)
    variance = (out["residual_sq_sum"] - (out["residual_sum"] ** 2) / n) / denom
    out["residual_std"] = np.sqrt(variance.clip(lower=0))
    return out


def add_site_reductions(overall: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, sub in overall.groupby("target"):
        base = sub[sub["model"].eq("mf2013_basic")]
        if base.empty:
            continue
        base_row = base.iloc[0]
        for _, row in sub.iterrows():
            item = row.to_dict()
            item["mae_reduction_vs_basic"] = base_row["mae"] - row["mae"]
            item["percent_mae_reduction_vs_basic"] = (
                (base_row["mae"] - row["mae"]) / base_row["mae"] * 100.0
            )
            rows.append(item)
    return pd.DataFrame(rows)


def station_site_correlations(station_terms: pd.DataFrame, min_records: int) -> pd.DataFrame:
    features = [
        ("log10(VS10)", "vs10"),
        ("log10(VS20)", "vs20"),
        ("log10(VS30)", "vs30"),
        ("log10(AVS30)", "avs30"),
        ("log10(D1100)", "d1100"),
        ("log10(D1400)", "d1400"),
        ("log10(Dbase)", "dbase"),
    ]
    rows = []
    supported = station_terms[station_terms["n_records"].ge(min_records)].copy()
    for (target, target_label, model), sub in supported.groupby(["target", "target_label", "model"]):
        for feature_label, col in features:
            x = pd.to_numeric(sub[col], errors="coerce")
            y = pd.to_numeric(sub["mean_residual"], errors="coerce")
            valid = x.gt(0) & y.notna()
            if valid.sum() < 30:
                continue
            rho, pvalue = stats.spearmanr(np.log10(x.loc[valid]), y.loc[valid])
            rows.append(
                {
                    "target": target,
                    "target_label": target_label,
                    "model": model,
                    "feature": feature_label,
                    "n_sites": int(valid.sum()),
                    "min_records_per_site": min_records,
                    "spearman_rho": float(rho),
                    "pvalue": float(pvalue),
                }
            )
    return pd.DataFrame(rows)


def site_correction_effects(overall: pd.DataFrame, correlations: pd.DataFrame) -> pd.DataFrame:
    metric_rows = []
    for _, row in add_site_reductions(overall).iterrows():
        if row["model"] != "mf2013_site":
            continue
        basic = overall[(overall["target"].eq(row["target"])) & overall["model"].eq("mf2013_basic")].iloc[0]
        metric_rows.append(
            {
                "target": row["target"],
                "target_label": row["target_label"],
                "effect_type": "mae",
                "feature": "overall",
                "basic_value": basic["mae"],
                "site_value": row["mae"],
                "absolute_reduction": basic["mae"] - row["mae"],
                "relative_reduction_pct": row["percent_mae_reduction_vs_basic"],
                "n": int(row["n_records"]),
            }
        )
    corr_rows = []
    if not correlations.empty:
        for (target, target_label, feature), sub in correlations.groupby(["target", "target_label", "feature"]):
            basic = sub[sub["model"].eq("mf2013_basic")]
            site = sub[sub["model"].eq("mf2013_site")]
            if basic.empty or site.empty:
                continue
            basic_rho = float(basic.iloc[0]["spearman_rho"])
            site_rho = float(site.iloc[0]["spearman_rho"])
            corr_rows.append(
                {
                    "target": target,
                    "target_label": target_label,
                    "effect_type": "station_residual_spearman_absrho",
                    "feature": feature,
                    "basic_value": basic_rho,
                    "site_value": site_rho,
                    "absolute_reduction": abs(basic_rho) - abs(site_rho),
                    "relative_reduction_pct": (
                        (abs(basic_rho) - abs(site_rho)) / abs(basic_rho) * 100.0
                        if abs(basic_rho) > 0
                        else np.nan
                    ),
                    "n": int(min(basic.iloc[0]["n_sites"], site.iloc[0]["n_sites"])),
                }
            )
    return pd.DataFrame(metric_rows + corr_rows)


def plot_outputs(overall: pd.DataFrame, correlations: pd.DataFrame, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    labels = [item["target_label"] for item in TARGETS]
    pivot = overall.pivot(index="target_label", columns="model", values="mae")
    pivot = pivot.loc[[label for label in labels if label in pivot.index]]
    x = np.arange(len(pivot.index))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.bar(x - width / 2, pivot["mf2013_basic"], width, label="MF2013 basic")
    ax.bar(x + width / 2, pivot["mf2013_site"], width, label="MF2013 + D1400/AVS30")
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=15, ha="right")
    ax.set_ylabel("MAE in log10 units")
    ax.set_title("Official J-SHIS Flatfile Residuals vs MF2013")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "jshis_mf2013_official_mae_by_target.png", dpi=180)
    plt.close(fig)

    if not correlations.empty:
        corr = correlations.copy()
        corr["column"] = corr["model"].str.replace("mf2013_", "", regex=False) + " " + corr["feature"]
        heat = corr.pivot(index="target_label", columns="column", values="spearman_rho")
        heat = heat.loc[[label for label in labels if label in heat.index]]
        vmax = np.nanmax(np.abs(heat.to_numpy()))
        fig, ax = plt.subplots(figsize=(12.5, 5.2))
        im = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels(heat.index)
        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels(heat.columns, rotation=35, ha="right", fontsize=8)
        for i in range(heat.shape[0]):
            for j in range(heat.shape[1]):
                val = heat.iloc[i, j]
                if np.isfinite(val):
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
        ax.set_title("Station Mean Residuals vs Official Site Parameters")
        fig.colorbar(im, ax=ax, label="Spearman rho")
        fig.tight_layout()
        fig.savefig(figure_dir / "jshis_mf2013_official_station_site_correlation_heatmap.png", dpi=180)
        plt.close(fig)


def write_report(
    output_dir: Path,
    zip_path: Path,
    total_smrec_rows: int,
    joined_rows: int,
    overall: pd.DataFrame,
    by_source: pd.DataFrame,
    station_terms: pd.DataFrame,
    correlations: pd.DataFrame,
    skipped_summary: pd.DataFrame,
    effects: pd.DataFrame,
    min_station_records: int,
) -> None:
    lines = ["# Official J-SHIS Flatfile Residuals Against MF2013 Core Equations\n\n"]
    lines.append("## Source data\n")
    lines.append(f"- Archive: `{zip_path.name}`.\n")
    lines.append("- Official flatfile schemas used: `smrec_schema.tsv`, `source_schema.tsv`, `site_schema.tsv`.\n")
    lines.append("- MF2013 coefficients: official J-SHIS/NIED `MF13rev_coefs.csv`.\n")
    lines.append("- Source URLs: `https://www.j-shis.bosai.go.jp/en/labs/ground-motion-flatfile/` and `https://www.j-shis.bosai.go.jp/en/labs/mf2013/`.\n\n")

    lines.append("## Coverage\n")
    lines.append(f"- `smrec_schema.tsv` rows read: {total_smrec_rows:,}.\n")
    lines.append(f"- Rows joined to source and site schema: {joined_rows:,}.\n")
    for _, row in skipped_summary.iterrows():
        lines.append(
            f"- {row['target_label']} / {row['model']}: usable residual rows={int(row['usable_rows']):,}, "
            f"missing/invalid rows={int(row['invalid_rows']):,}.\n"
        )
    lines.append("\n## Overall residual metrics\n")
    overall_reduced = add_site_reductions(overall)
    for _, row in overall_reduced.sort_values(["target_label", "model"]).iterrows():
        lines.append(
            f"- {row['target_label']} / {row['model']}: n={int(row['n_records']):,}, "
            f"MAE={row['mae']:.4f}, RMSE={row['rmse']:.4f}, "
            f"bias observed-minus-predicted={row['mean_residual']:.4f}, "
            f"MAE reduction vs basic={row['percent_mae_reduction_vs_basic']:.1f}%.\n"
        )
    lines.append("\n## Source-class coverage\n")
    source_pivot = by_source.pivot_table(
        index=["target_label", "model"],
        columns="source_class_label",
        values="n_records",
        aggfunc="sum",
        fill_value=0,
    )
    for idx, row in source_pivot.iterrows():
        parts = ", ".join(f"{col}={int(val):,}" for col, val in row.items())
        lines.append(f"- {idx[0]} / {idx[1]}: {parts}.\n")

    lines.append("\n## Strongest station-residual site associations\n")
    if correlations.empty:
        lines.append("- No station-site correlations met the support threshold.\n")
    else:
        best = correlations.reindex(correlations["spearman_rho"].abs().sort_values(ascending=False).index).head(12)
        for _, row in best.iterrows():
            lines.append(
                f"- {row['target_label']} / {row['model']} vs {row['feature']}: "
                f"rho={row['spearman_rho']:.3f}, n_sites={int(row['n_sites'])}.\n"
            )

    lines.append("\n## Site-correction effect summary\n")
    mae_effects = effects[effects["effect_type"].eq("mae")]
    for _, row in mae_effects.iterrows():
        lines.append(
            f"- {row['target_label']}: MF2013 site terms reduce MAE by "
            f"{row['relative_reduction_pct']:.1f}% ({row['basic_value']:.4f} to {row['site_value']:.4f}).\n"
        )
    corr_effects = effects[effects["effect_type"].eq("station_residual_spearman_absrho")]
    if not corr_effects.empty:
        best_corr_effects = corr_effects.sort_values("absolute_reduction", ascending=False).head(8)
        for _, row in best_corr_effects.iterrows():
            lines.append(
                f"- {row['target_label']} / {row['feature']}: |rho| reduction "
                f"{row['absolute_reduction']:.3f} ({row['relative_reduction_pct']:.1f}%), "
                f"basic rho={row['basic_value']:.3f}, site rho={row['site_value']:.3f}.\n"
            )

    supported_sites = station_terms[station_terms["n_records"].ge(min_station_records)]
    lines.append("\n## Station residual support\n")
    lines.append(
        f"- Station terms with at least {min_station_records} records per target/model: {len(supported_sites):,} rows across "
        f"{supported_sites['siteid2'].nunique():,} site IDs.\n"
    )

    lines.append("\n## Outputs\n")
    for name in [
        "jshis_mf2013_official_overall_metrics.csv",
        "jshis_mf2013_official_by_network_metrics.csv",
        "jshis_mf2013_official_by_source_class_metrics.csv",
        "jshis_mf2013_official_station_residual_terms.csv",
        "jshis_mf2013_official_event_residual_terms.csv",
        "jshis_mf2013_official_station_site_correlations.csv",
        "jshis_mf2013_official_site_correction_effects.csv",
        "jshis_mf2013_official_residual_sample.csv",
        "figures/jshis_mf2013_official_mae_by_target.png",
        "figures/jshis_mf2013_official_station_site_correlation_heatmap.png",
    ]:
        lines.append(f"- `{name}`\n")

    lines.append("\n## Guardrails\n")
    lines.append("- This is an official-data MF2013 core-equation residual analysis, not a full J-SHIS hazard-map reproduction.\n")
    lines.append("- The `AI` anomalous-intensity and `PH` Philippine Sea Plate correction terms are not applied here.\n")
    lines.append("- Residuals are raw observed-minus-predicted log10 values, not train-calibrated predictions.\n")
    lines.append("- `flatfile_sub1-v2024.zip` is the M>=5 / fault-distance<=300 km public subset, not the full 2.0 GB flatfile.\n")
    (output_dir / "jshis_mf2013_official_residuals.md").write_text("".join(lines), encoding="utf-8")


def run_analysis(
    zip_path: Path,
    coeff_path: Path,
    output_dir: Path,
    chunksize: int,
    sample_rows_per_chunk: int,
    min_station_records: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    coeffs = load_coeffs(coeff_path)
    site, source = read_schema_tables(zip_path)

    overall_parts: list[pd.DataFrame] = []
    network_parts: list[pd.DataFrame] = []
    source_parts: list[pd.DataFrame] = []
    station_parts: list[pd.DataFrame] = []
    event_parts: list[pd.DataFrame] = []
    residual_samples: list[pd.DataFrame] = []
    skipped_rows: list[dict[str, object]] = []

    total_smrec_rows = 0
    joined_rows = 0
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
                total_smrec_rows += len(chunk)
                chunk["site_id"] = pd.to_numeric(chunk["site_id"], errors="coerce").astype("Int64")
                chunk["siteid2"] = (chunk["site_id"] // 10).astype("Int64")
                chunk["eq_source_id"] = pd.to_numeric(chunk["eq_source_id"], errors="coerce").astype("Int64")
                joined = chunk.merge(site, on="siteid2", how="left").merge(source, on="eq_source_id", how="left")
                joined_rows += len(joined)

                residual_frames = []
                for target_info in TARGETS:
                    for model_name, include_site in [
                        ("mf2013_basic", False),
                        ("mf2013_site", True),
                    ]:
                        residuals = residual_frame(joined, coeffs, target_info, model_name, include_site)
                        residual_frames.append(residuals)
                        skipped_rows.append(
                            {
                                "chunk": chunk_index,
                                "target": target_info["target"],
                                "target_label": target_info["target_label"],
                                "model": model_name,
                                "input_rows": len(joined),
                                "usable_rows": len(residuals),
                                "invalid_rows": len(joined) - len(residuals),
                            }
                        )
                residual_long = pd.concat([frame for frame in residual_frames if not frame.empty], ignore_index=True)
                residual_long["source_class_label"] = residual_long["eq_location_type_id"].map(SOURCE_CLASS_LABELS)
                residual_long["network_label"] = residual_long["obs_network_id"].map(NETWORK_LABELS)

                overall_parts.append(aggregate_residuals(residual_long, ["target", "target_label", "period", "model"]))
                network_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model", "obs_network_id", "network_label"],
                    )
                )
                source_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model", "eq_location_type_id", "source_class_label"],
                    )
                )
                station_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model", "siteid2"],
                    )
                )
                event_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model", "eq_source_id"],
                    )
                )

                if sample_rows_per_chunk > 0 and not residual_long.empty:
                    sample = residual_long.sample(
                        n=min(sample_rows_per_chunk, len(residual_long)),
                        random_state=20260605 + chunk_index,
                    )
                    residual_samples.append(sample)

                print(
                    f"chunk={chunk_index} smrec_rows={total_smrec_rows:,} "
                    f"residual_rows_this_chunk={len(residual_long):,}",
                    flush=True,
                )

    overall = combine_aggs(overall_parts, ["target", "target_label", "period", "model"])
    by_network = combine_aggs(
        network_parts,
        ["target", "target_label", "period", "model", "obs_network_id", "network_label"],
    )
    by_source = combine_aggs(
        source_parts,
        ["target", "target_label", "period", "model", "eq_location_type_id", "source_class_label"],
    )
    station_terms = combine_aggs(station_parts, ["target", "target_label", "period", "model", "siteid2"])
    event_terms = combine_aggs(event_parts, ["target", "target_label", "period", "model", "eq_source_id"])

    station_terms = station_terms.merge(site, on="siteid2", how="left")
    event_terms = event_terms.merge(source, on="eq_source_id", how="left")
    correlations = station_site_correlations(station_terms, min_records=min_station_records)
    effects = site_correction_effects(overall, correlations)

    skipped = pd.DataFrame(skipped_rows)
    skipped_summary = (
        skipped.groupby(["target", "target_label", "model"], as_index=False)[["usable_rows", "invalid_rows"]]
        .sum()
        .sort_values(["target_label", "model"])
    )

    overall.to_csv(output_dir / "jshis_mf2013_official_overall_metrics.csv", index=False)
    by_network.to_csv(output_dir / "jshis_mf2013_official_by_network_metrics.csv", index=False)
    by_source.to_csv(output_dir / "jshis_mf2013_official_by_source_class_metrics.csv", index=False)
    station_terms.to_csv(output_dir / "jshis_mf2013_official_station_residual_terms.csv", index=False)
    event_terms.to_csv(output_dir / "jshis_mf2013_official_event_residual_terms.csv", index=False)
    correlations.to_csv(output_dir / "jshis_mf2013_official_station_site_correlations.csv", index=False)
    effects.to_csv(output_dir / "jshis_mf2013_official_site_correction_effects.csv", index=False)
    skipped_summary.to_csv(output_dir / "jshis_mf2013_official_coverage_summary.csv", index=False)
    if residual_samples:
        pd.concat(residual_samples, ignore_index=True).to_csv(
            output_dir / "jshis_mf2013_official_residual_sample.csv",
            index=False,
        )

    plot_outputs(overall, correlations, output_dir)
    write_report(
        output_dir=output_dir,
        zip_path=zip_path,
        total_smrec_rows=total_smrec_rows,
        joined_rows=joined_rows,
        overall=overall,
        by_source=by_source,
        station_terms=station_terms,
        correlations=correlations,
        skipped_summary=skipped_summary,
        effects=effects,
        min_station_records=min_station_records,
    )

    print(f"total_smrec_rows={total_smrec_rows:,}")
    print(f"overall_metric_rows={len(overall):,}")
    print(f"station_term_rows={len(station_terms):,}")
    print(f"event_term_rows={len(event_terms):,}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--coeffs", type=Path, default=DEFAULT_COEFFS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--sample-rows-per-chunk", type=int, default=5_000)
    parser.add_argument("--min-station-records", type=int, default=20)
    args = parser.parse_args()
    run_analysis(
        zip_path=args.zip,
        coeff_path=args.coeffs,
        output_dir=args.output_dir,
        chunksize=args.chunksize,
        sample_rows_per_chunk=args.sample_rows_per_chunk,
        min_station_records=args.min_station_records,
    )


if __name__ == "__main__":
    main()
