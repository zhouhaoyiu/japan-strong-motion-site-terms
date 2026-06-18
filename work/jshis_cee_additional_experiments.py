#!/usr/bin/env python3
"""Additional CEE submission checks for the J-SHIS station-residual study.

The script writes audit tables that support the manuscript and supplementary
files. It starts from existing derived station-level products and avoids
re-reading the full waveform/record flatfile.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from jshis_nonergodic_station_residual_model import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_SITE_SCHEMA,
    DEFAULT_STATION_TERMS,
    TARGET_ORDER,
    load_station_dataset,
    make_spatial_folds,
    weighted_mae,
    weighted_rmse,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TODAY = "2026-06-14"

FEATURE_SETS: dict[str, list[str]] = {
    "mf2013_d1400_avs30": ["log_d1400", "log_avs30"],
    "d1400_only": ["log_d1400"],
    "avs30_only": ["log_avs30"],
    "vs20_dbase": ["log_vs20", "log_dbase"],
    "official_four": ["log_d1400", "log_avs30", "log_vs20", "log_dbase"],
    "depth_family": ["log_d1100", "log_d1400", "log_d1700", "log_d2100", "log_dbase"],
    "site_no_location": [
        "log_vs10",
        "log_vs20",
        "log_vs30",
        "log_avs30",
        "log_d1100",
        "log_d1400",
        "log_d1700",
        "log_d2100",
        "log_dbase",
        "elevation",
        "sensor_depth_glminus",
        "dist_vf_mf13_nejapan",
        "dist_vf_mf13_swjapan",
    ],
    "full_site_location": [
        "log_vs10",
        "log_vs20",
        "log_vs30",
        "log_avs30",
        "log_d1100",
        "log_d1400",
        "log_d1700",
        "log_d2100",
        "log_dbase",
        "elevation",
        "sensor_depth_glminus",
        "dist_vf_mf13_nejapan",
        "dist_vf_mf13_swjapan",
        "lon",
        "lat",
    ],
}

CATEGORICAL_SETS: dict[str, list[str]] = {
    "mf2013_d1400_avs30": [],
    "d1400_only": [],
    "avs30_only": [],
    "vs20_dbase": [],
    "official_four": [],
    "depth_family": [],
    "site_no_location": ["network_label", "installation_situation_id"],
    "full_site_location": ["network_label", "installation_situation_id"],
}

REGION_PREFIXES: dict[str, str] = {
    # Hokkaido subprefecture and NIED-style prefixes.
    "ABS": "Hokkaido",
    "HDK": "Hokkaido",
    "HKD": "Hokkaido",
    "HYM": "Hokkaido",
    "IBU": "Hokkaido",
    "IKR": "Hokkaido",
    "KKW": "Hokkaido",
    "KSR": "Hokkaido",
    "NMR": "Hokkaido",
    "OSM": "Hokkaido",
    "RMI": "Hokkaido",
    "SBS": "Hokkaido",
    "SOY": "Hokkaido",
    "SRC": "Hokkaido",
    "TKC": "Hokkaido",
    # Tohoku.
    "AKT": "Tohoku",
    "AOM": "Tohoku",
    "FKS": "Tohoku",
    "IWT": "Tohoku",
    "MYG": "Tohoku",
    "YMT": "Tohoku",
    # Kanto.
    "CHB": "Kanto",
    "GNM": "Kanto",
    "IBR": "Kanto",
    "KNG": "Kanto",
    "SIT": "Kanto",
    "TCG": "Kanto",
    "TKY": "Kanto",
    # Chubu.
    "AIC": "Chubu",
    "FKI": "Chubu",
    "GIF": "Chubu",
    "ISK": "Chubu",
    "NGN": "Chubu",
    "NIG": "Chubu",
    "SZO": "Chubu",
    "TYM": "Chubu",
    "YMN": "Chubu",
    # Kinki.
    "HYG": "Kinki",
    "KYT": "Kinki",
    "MIE": "Kinki",
    "NAR": "Kinki",
    "OSK": "Kinki",
    "SIG": "Kinki",
    "WKY": "Kinki",
    # Chugoku and Shikoku.
    "EHM": "Chugoku_Shikoku",
    "HRS": "Chugoku_Shikoku",
    "KGW": "Chugoku_Shikoku",
    "KOC": "Chugoku_Shikoku",
    "OKY": "Chugoku_Shikoku",
    "SMN": "Chugoku_Shikoku",
    "TKS": "Chugoku_Shikoku",
    "TTR": "Chugoku_Shikoku",
    "YMG": "Chugoku_Shikoku",
    # Kyushu and Okinawa.
    "FKO": "Kyushu_Okinawa",
    "KGS": "Kyushu_Okinawa",
    "KMM": "Kyushu_Okinawa",
    "MYZ": "Kyushu_Okinawa",
    "NGS": "Kyushu_Okinawa",
    "OIT": "Kyushu_Okinawa",
    "OKN": "Kyushu_Okinawa",
    "SAG": "Kyushu_Okinawa",
}

TARGET_SHORT = {
    "PGA RotD50": "PGA",
    "SA(0.3s) RotD50": "SA(0.3 s)",
    "SA(1.0s) RotD50": "SA(1.0 s)",
    "SA(3.0s) RotD50": "SA(3.0 s)",
}


def positive_log10(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return np.log10(values.where(values > 0))


def format_md_cell(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (float, np.floating)):
        if math.isinf(float(value)):
            return "inf"
        return f"{float(value):.4g}"
    text = str(value).replace("\n", " ").replace("|", "\\|")
    return text


def df_to_markdown(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows._"
    columns = [str(col) for col in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(format_md_cell(row[col]) for col in frame.columns) + " |")
    return "\n".join(lines)


def station_prefix(site_code: object) -> str:
    text = str(site_code).strip().upper()
    match = re.match(r"([A-Z]+)", text)
    if not match:
        return ""
    return match.group(1)[:3]


def assign_region(site_code: object, lon: object, lat: object) -> str:
    prefix = station_prefix(site_code)
    if prefix in REGION_PREFIXES:
        return REGION_PREFIXES[prefix]
    lon_value = pd.to_numeric(pd.Series([lon]), errors="coerce").iloc[0]
    lat_value = pd.to_numeric(pd.Series([lat]), errors="coerce").iloc[0]
    if pd.isna(lon_value) or pd.isna(lat_value):
        return "Unknown"
    if lat_value >= 41.0:
        return "Hokkaido"
    if lat_value >= 37.0 and lon_value >= 139.0:
        return "Tohoku"
    if 138.0 <= lon_value <= 141.0 and 34.0 <= lat_value < 37.2:
        return "Kanto"
    if 136.0 <= lon_value < 139.5 and 34.0 <= lat_value < 38.0:
        return "Chubu"
    if 134.0 <= lon_value < 136.8 and 33.0 <= lat_value < 36.2:
        return "Kinki"
    if 131.0 <= lon_value < 134.8 and 32.0 <= lat_value < 36.0:
        return "Chugoku_Shikoku"
    return "Kyushu_Okinawa"


def add_region(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["station_prefix"] = out["site_code"].map(station_prefix)
    out["region"] = [
        assign_region(site_code, lon, lat)
        for site_code, lon, lat in zip(out["site_code"], out["lon"], out["lat"], strict=True)
    ]
    return out


def estimator(seed: int) -> HistGradientBoostingRegressor:
    return HistGradientBoostingRegressor(
        loss="squared_error",
        learning_rate=0.045,
        max_iter=180,
        max_leaf_nodes=13,
        min_samples_leaf=35,
        l2_regularization=0.1,
        random_state=seed,
    )


def build_pipeline(feature_set: str, seed: int) -> Pipeline:
    numeric = FEATURE_SETS[feature_set]
    categorical = CATEGORICAL_SETS[feature_set]
    transformers: list[tuple[str, Pipeline, list[str]]] = [
        ("num", Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric)
    ]
    if categorical:
        transformers.append(
            (
                "cat",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                    ]
                ),
                categorical,
            )
        )
    return Pipeline([("preprocess", ColumnTransformer(transformers)), ("model", estimator(seed))])


def metric_row(
    target: str,
    split: str,
    fold: object,
    model: str,
    train_n: int,
    test_n: int,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    weights: np.ndarray,
    baseline_rmse: float,
    baseline_mae: float,
) -> dict:
    mae = weighted_mae(y_true, y_pred, weights)
    rmse = weighted_rmse(y_true, y_pred, weights)
    return {
        "target_label": target,
        "period_label": TARGET_SHORT[target],
        "split": split,
        "fold": fold,
        "model": model,
        "n_train_sites": int(train_n),
        "n_test_sites": int(test_n),
        "test_record_weight": float(weights.sum()),
        "weighted_mae": mae,
        "weighted_rmse": rmse,
        "mae_reduction_pct": 100.0 * (baseline_mae - mae) / baseline_mae if baseline_mae else math.nan,
        "rmse_reduction_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse if baseline_rmse else math.nan,
        "weighted_r2_vs_zero": 1.0 - (rmse**2 / baseline_rmse**2) if baseline_rmse else math.nan,
    }


def summarize_metrics(fold_metrics: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    summary = (
        fold_metrics.groupby(group_cols, observed=True)
        .agg(
            folds=("fold", "nunique"),
            n_test_sites_total=("n_test_sites", "sum"),
            n_test_sites_mean=("n_test_sites", "mean"),
            weighted_mae_mean=("weighted_mae", "mean"),
            weighted_mae_sd=("weighted_mae", "std"),
            weighted_rmse_mean=("weighted_rmse", "mean"),
            weighted_rmse_sd=("weighted_rmse", "std"),
            mae_reduction_pct_mean=("mae_reduction_pct", "mean"),
            mae_reduction_pct_min=("mae_reduction_pct", "min"),
            mae_reduction_pct_max=("mae_reduction_pct", "max"),
            rmse_reduction_pct_mean=("rmse_reduction_pct", "mean"),
            rmse_reduction_pct_min=("rmse_reduction_pct", "min"),
            rmse_reduction_pct_max=("rmse_reduction_pct", "max"),
            weighted_r2_vs_zero_mean=("weighted_r2_vs_zero", "mean"),
        )
        .reset_index()
    )
    summary["target_label"] = pd.Categorical(summary["target_label"], TARGET_ORDER, ordered=True)
    return summary.sort_values(group_cols).reset_index(drop=True)


def run_region_leaveout(frame: pd.DataFrame, output_dir: Path, seed: int) -> None:
    rows: list[dict] = []
    models = ["mf2013_d1400_avs30", "official_four", "site_no_location", "full_site_location"]
    for target in TARGET_ORDER:
        target_frame = frame[frame["target_label"].astype(str).eq(target)].copy().reset_index(drop=True)
        y = target_frame["mean_residual"].to_numpy(dtype=float)
        weights = target_frame["n_records"].to_numpy(dtype=float)
        for region in sorted(target_frame["region"].dropna().unique()):
            test_idx = np.flatnonzero(target_frame["region"].to_numpy() == region)
            train_idx = np.flatnonzero(target_frame["region"].to_numpy() != region)
            if len(test_idx) < 20 or len(train_idx) < 100:
                continue
            y_test = y[test_idx]
            w_test = weights[test_idx]
            baseline = np.zeros_like(y_test)
            baseline_mae = weighted_mae(y_test, baseline, w_test)
            baseline_rmse = weighted_rmse(y_test, baseline, w_test)
            rows.append(
                metric_row(
                    target,
                    "region_leaveout",
                    region,
                    "zero_station_residual",
                    len(train_idx),
                    len(test_idx),
                    y_test,
                    baseline,
                    w_test,
                    baseline_rmse,
                    baseline_mae,
                )
            )
            for feature_set in models:
                pipe = build_pipeline(feature_set, seed)
                pipe.fit(target_frame.iloc[train_idx], y[train_idx], model__sample_weight=weights[train_idx])
                pred = pipe.predict(target_frame.iloc[test_idx])
                rows.append(
                    metric_row(
                        target,
                        "region_leaveout",
                        region,
                        feature_set,
                        len(train_idx),
                        len(test_idx),
                        y_test,
                        pred,
                        w_test,
                        baseline_rmse,
                        baseline_mae,
                    )
                )

    fold_metrics = pd.DataFrame(rows)
    fold_metrics.to_csv(output_dir / "jshis_station_region_leaveout_fold_metrics.csv", index=False)
    summary = summarize_metrics(fold_metrics, ["target_label", "period_label", "model", "split"])
    summary.to_csv(output_dir / "jshis_station_region_leaveout_summary.csv", index=False)


def run_variable_substitution(frame: pd.DataFrame, output_dir: Path, seed: int) -> None:
    rows: list[dict] = []
    for target in TARGET_ORDER:
        target_frame = frame[frame["target_label"].astype(str).eq(target)].copy().reset_index(drop=True)
        folds = make_spatial_folds(target_frame, n_splits=5, seed=seed)
        y = target_frame["mean_residual"].to_numpy(dtype=float)
        weights = target_frame["n_records"].to_numpy(dtype=float)
        for fold in sorted(np.unique(folds)):
            test_idx = np.flatnonzero(folds == fold)
            train_idx = np.flatnonzero(folds != fold)
            y_test = y[test_idx]
            w_test = weights[test_idx]
            baseline = np.zeros_like(y_test)
            baseline_mae = weighted_mae(y_test, baseline, w_test)
            baseline_rmse = weighted_rmse(y_test, baseline, w_test)
            rows.append(
                metric_row(
                    target,
                    "spatial_block",
                    int(fold),
                    "zero_station_residual",
                    len(train_idx),
                    len(test_idx),
                    y_test,
                    baseline,
                    w_test,
                    baseline_rmse,
                    baseline_mae,
                )
            )
            for feature_set in FEATURE_SETS:
                pipe = build_pipeline(feature_set, seed)
                pipe.fit(target_frame.iloc[train_idx], y[train_idx], model__sample_weight=weights[train_idx])
                pred = pipe.predict(target_frame.iloc[test_idx])
                rows.append(
                    metric_row(
                        target,
                        "spatial_block",
                        int(fold),
                        feature_set,
                        len(train_idx),
                        len(test_idx),
                        y_test,
                        pred,
                        w_test,
                        baseline_rmse,
                        baseline_mae,
                    )
                )

    fold_metrics = pd.DataFrame(rows)
    fold_metrics.to_csv(output_dir / "jshis_station_variable_substitution_fold_metrics.csv", index=False)
    summary = summarize_metrics(fold_metrics, ["target_label", "period_label", "model", "split"])
    summary.to_csv(output_dir / "jshis_station_variable_substitution_summary.csv", index=False)


def run_collinearity(frame: pd.DataFrame, output_dir: Path) -> None:
    site_frame = frame[frame["target_label"].astype(str).eq("SA(3.0s) RotD50")].drop_duplicates("siteid2").copy()
    variables = ["log_vs20", "log_vs30", "log_avs30", "log_d1100", "log_d1400", "log_d2100", "log_dbase"]
    labels = {
        "log_vs20": "VS20",
        "log_vs30": "VS30",
        "log_avs30": "AVS30",
        "log_d1100": "D1100",
        "log_d1400": "D1400",
        "log_d2100": "D2100",
        "log_dbase": "Dbase",
    }
    corr = site_frame[variables].corr(method="spearman")
    rows = []
    for i, left in enumerate(variables):
        for right in variables[i + 1 :]:
            rows.append(
                {
                    "left_variable": labels[left],
                    "right_variable": labels[right],
                    "spearman_rho": float(corr.loc[left, right]),
                    "n_complete": int(site_frame[[left, right]].dropna().shape[0]),
                }
            )
    pd.DataFrame(rows).to_csv(output_dir / "jshis_station_site_variable_collinearity_pairs.csv", index=False)

    vif_vars = ["log_vs20", "log_avs30", "log_d1400", "log_dbase"]
    matrix = site_frame[vif_vars].copy()
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.fillna(matrix.median(numeric_only=True))
    vif_rows = []
    for variable in vif_vars:
        others = [col for col in vif_vars if col != variable]
        model = LinearRegression().fit(matrix[others], matrix[variable])
        r2 = float(model.score(matrix[others], matrix[variable]))
        vif = math.inf if r2 >= 1.0 else 1.0 / (1.0 - r2)
        vif_rows.append({"variable": labels[variable], "r2_from_other_variables": r2, "vif": vif})
    pd.DataFrame(vif_rows).to_csv(output_dir / "jshis_station_site_variable_vif.csv", index=False)


def run_positive_correction_audit(frame: pd.DataFrame, output_dir: Path) -> None:
    values_path = output_dir / "jshis_official_response_sa3_station_values.csv"
    if not values_path.exists():
        return
    values = pd.read_csv(values_path)
    sa3 = frame[frame["target_label"].astype(str).eq("SA(3.0s) RotD50")][
        ["siteid2", "site_code", "network_label", "n_records", "mean_residual", "region", "lon", "lat"]
    ].drop_duplicates("siteid2")
    values = values[values["probability_level"].eq("50y_10pct")].copy()
    values = values.merge(sa3, on=["siteid2", "site_code", "network_label", "n_records"], how="left", validate="one_to_one")
    values["positive_correction"] = values["station_multiplier"].gt(1.0)
    values["negative_correction"] = values["station_multiplier"].lt(1.0)
    values["abs_log10_correction"] = values["sa3_correction_log10"].abs()
    values = values.sort_values("station_multiplier", ascending=False)
    values.to_csv(output_dir / "jshis_positive_sa3_station_correction_audit.csv", index=False)

    positive = values[values["positive_correction"]].copy()
    summary_rows = [
        {
            "metric": "all_stations",
            "value": int(values.shape[0]),
            "unit": "stations",
            "note": "Stations with held-out SA(3.0 s) correction and official response-map match.",
        },
        {
            "metric": "positive_correction_stations",
            "value": int(positive.shape[0]),
            "unit": "stations",
            "note": "Station multiplier greater than 1.0 at 50-year 10% probability.",
        },
        {
            "metric": "positive_correction_share_pct",
            "value": 100.0 * positive.shape[0] / values.shape[0] if len(values) else math.nan,
            "unit": "percent",
            "note": "Share of matched stations with upward correction.",
        },
        {
            "metric": "positive_multiplier_median",
            "value": float(positive["station_multiplier"].median()) if len(positive) else math.nan,
            "unit": "ratio",
            "note": "Median upward multiplier among positive-correction stations.",
        },
        {
            "metric": "positive_multiplier_95pct",
            "value": float(positive["station_multiplier"].quantile(0.95)) if len(positive) else math.nan,
            "unit": "ratio",
            "note": "95th percentile upward multiplier among positive-correction stations.",
        },
        {
            "metric": "positive_multiplier_max",
            "value": float(positive["station_multiplier"].max()) if len(positive) else math.nan,
            "unit": "ratio",
            "note": "Largest upward multiplier among positive-correction stations.",
        },
    ]
    region_counts = (
        values.groupby(["region", "positive_correction"], dropna=False)
        .size()
        .reset_index(name="n_stations")
        .sort_values(["region", "positive_correction"])
    )
    for region, group in region_counts.groupby("region", dropna=False):
        total = int(group["n_stations"].sum())
        pos = int(group.loc[group["positive_correction"].eq(True), "n_stations"].sum())
        summary_rows.append(
            {
                "metric": f"positive_share_{region}",
                "value": 100.0 * pos / total if total else math.nan,
                "unit": "percent",
                "note": f"{pos}/{total} matched stations in {region}.",
            }
        )
    pd.DataFrame(summary_rows).to_csv(output_dir / "jshis_positive_sa3_station_correction_summary.csv", index=False)

    top = values.head(12)[
        [
            "site_code",
            "region",
            "network_label",
            "n_records",
            "station_multiplier",
            "sa3_correction_log10",
            "official_sa3_g",
            "corrected_sa3_g",
            "mean_residual",
        ]
    ]
    lines = [
        "# Positive SA(3.0 s) station-correction audit",
        "",
        f"Generated on {TODAY}. The audit uses the official 50-year 10% response-map match and held-out SA(3.0 s) station corrections.",
        "",
        "## Main counts",
        "",
    ]
    for row in summary_rows[:6]:
        lines.append(f"- {row['metric']}: {row['value']:.4g} {row['unit']}. {row['note']}")
    lines.extend(["", "## Largest upward corrections", "", df_to_markdown(top), ""])
    (output_dir / "jshis_positive_sa3_station_correction_audit.md").write_text("\n".join(lines), encoding="utf-8")


def parse_summary_value(text: str, label: str) -> int | None:
    match = re.search(re.escape(label) + r"\s*:\s*([\d,]+)", text)
    if not match:
        return None
    return int(match.group(1).replace(",", ""))


def run_matching_missingness_audit(output_dir: Path, site_schema_path: Path, station_terms_path: Path) -> None:
    terms = pd.read_csv(station_terms_path)
    site = pd.read_csv(site_schema_path)
    summary_path = output_dir / "jshis_smrec_schema_sub1_summary.md"
    summary_text = summary_path.read_text(encoding="utf-8") if summary_path.exists() else ""
    total_records = parse_summary_value(summary_text, "Strong-motion records")
    unique_record_sites = parse_summary_value(summary_text, "Unique `siteid2` values after `site_id // 10` mapping")
    unmatched_site_ids = parse_summary_value(summary_text, "Unmatched `site_id` values observed")
    if total_records is None:
        coverage_path = output_dir / "jshis_mf2013_official_coverage_summary.csv"
        if coverage_path.exists():
            coverage = pd.read_csv(coverage_path)
            first = coverage.iloc[0]
            total_records = int(first["usable_rows"] + first["invalid_rows"])
    rows: list[dict] = []

    def add(check: str, count: float, denominator: float | None, source: str, note: str) -> None:
        rows.append(
            {
                "check": check,
                "count": count,
                "denominator": denominator,
                "percent": 100.0 * count / denominator if denominator else math.nan,
                "source": source,
                "note": note,
            }
        )

    if total_records is not None:
        add("strong_motion_records", total_records, None, "jshis_smrec_schema_sub1_summary.md", "Public flatfile records in sub1-v2024.")
    if unique_record_sites is not None:
        add(
            "unique_record_sites",
            unique_record_sites,
            None,
            "jshis_smrec_schema_sub1_summary.md",
            "Unique site_id//10 values in the strong-motion records.",
        )
    add("site_schema_rows", site.shape[0], None, "jshis_site_schema_v2024_sub1.csv", "Rows in the public site schema.")
    if total_records is not None:
        add(
            "record_site_schema_matched_rows",
            total_records,
            total_records,
            "jshis_smrec_schema_sub1_summary.md",
            "Record-to-site match reported by the schema extraction audit.",
        )
    if unmatched_site_ids is not None:
        add(
            "unmatched_site_id_values",
            unmatched_site_ids,
            unique_record_sites,
            "jshis_smrec_schema_sub1_summary.md",
            "Site identifiers in records without a site-schema match.",
        )

    for variable in ["vs10", "vs20", "vs30", "avs30", "d1100", "d1400", "d1700", "d2100", "dbase"]:
        values = pd.to_numeric(site[variable], errors="coerce")
        valid = int(values.gt(0).sum())
        add(
            f"valid_{variable}",
            valid,
            int(site.shape[0]),
            "jshis_site_schema_v2024_sub1.csv",
            f"Positive finite {variable.upper()} values in the public site schema.",
        )

    coverage_path = output_dir / "jshis_mf2013_official_coverage_summary.csv"
    if coverage_path.exists():
        coverage = pd.read_csv(coverage_path)
        for _, row in coverage.iterrows():
            denominator = int(row["usable_rows"] + row["invalid_rows"])
            add(
                f"{row['target_label']}_{row['model']}_usable_records",
                int(row["usable_rows"]),
                denominator,
                "jshis_mf2013_official_coverage_summary.csv",
                "Rows with valid observed and predicted amplitudes for station residual aggregation.",
            )

    for (target, model), group in terms.groupby(["target_label", "model"], observed=True):
        add(
            f"{target}_{model}_station_terms",
            int(group.shape[0]),
            unique_record_sites,
            "jshis_mf2013_official_station_residual_terms.csv",
            "Station-mean residual terms before the minimum-record threshold.",
        )
        enough = int(group["n_records"].ge(20).sum())
        add(
            f"{target}_{model}_station_terms_min20",
            enough,
            int(group.shape[0]),
            "jshis_mf2013_official_station_residual_terms.csv",
            "Station residual terms retained for held-out station modeling.",
        )

    response_values = output_dir / "jshis_official_response_sa3_station_values.csv"
    if response_values.exists():
        values = pd.read_csv(response_values)
        for probability, group in values.groupby("probability_level"):
            add(
                f"official_sa3_response_match_{probability}",
                int(group.shape[0]),
                int(values["siteid2"].nunique()),
                "jshis_official_response_sa3_station_values.csv",
                "Held-out station sites matched to official 250 m response-spectrum map values.",
            )

    audit = pd.DataFrame(rows)
    audit.to_csv(output_dir / "jshis_station_matching_missingness_audit.csv", index=False)

    key_rows = audit[audit["check"].isin(["strong_motion_records", "unique_record_sites", "site_schema_rows", "record_site_schema_matched_rows"])]
    lines = [
        "# Station matching and missingness audit",
        "",
        f"Generated on {TODAY}. Counts come from derived public-data audit tables.",
        "",
        df_to_markdown(key_rows),
        "",
        "The CSV file contains variable-level missingness, MF2013 residual coverage, minimum-record thresholds, and official SA(3.0 s) response-map matching.",
        "",
    ]
    (output_dir / "jshis_station_matching_missingness_audit.md").write_text("\n".join(lines), encoding="utf-8")


def run_provenance_release_manifest(output_dir: Path) -> None:
    rows = [
        {
            "asset": "J-SHIS/NIED strong-motion flatfile",
            "role": "Input ground-motion records",
            "version_or_date": "sub1-v2024",
            "source": "J-SHIS public download service, National Research Institute for Earth Science and Disaster Resilience",
            "source_url_or_identifier": "https://www.j-shis.bosai.go.jp/",
            "local_reference": "work/external_data/jshis_gmf/flatfile_sub1-v2024.zip",
            "derived_outputs": "jshis_smrec_schema_sub1_summary.md; jshis_mf2013_official_station_residual_terms.csv",
            "access_date": TODAY,
            "redistribution_note": "Record-level redistribution follows the original J-SHIS/NIED terms; submit derived audit tables and code separately.",
        },
        {
            "asset": "J-SHIS site schema",
            "role": "Official station/site parameters",
            "version_or_date": "sub1-v2024",
            "source": "site_schema.tsv in public flatfile package",
            "source_url_or_identifier": "https://www.j-shis.bosai.go.jp/",
            "local_reference": "outputs/jshis_site_schema_v2024_sub1.csv",
            "derived_outputs": "station matching, site-variable models, collinearity tables",
            "access_date": TODAY,
            "redistribution_note": "Use derived station-level summaries in the manuscript package.",
        },
        {
            "asset": "J-SHIS source schema",
            "role": "Earthquake source metadata",
            "version_or_date": "sub1-v2024",
            "source": "source_schema.tsv in public flatfile package",
            "source_url_or_identifier": "https://www.j-shis.bosai.go.jp/",
            "local_reference": "outputs/jshis_source_schema_v2024_sub1.csv",
            "derived_outputs": "event-exclusion and source-class sensitivity tables",
            "access_date": TODAY,
            "redistribution_note": "Use derived source-level summaries in the manuscript package.",
        },
        {
            "asset": "Morikawa and Fujiwara (2013) model implementation",
            "role": "Published ergodic reference model",
            "version_or_date": "2013 paper coefficients as implemented in project scripts",
            "source": "Morikawa and Fujiwara (2013), Bulletin of the Seismological Society of America",
            "source_url_or_identifier": "https://doi.org/10.1785/0120120097",
            "local_reference": "work/jshis_mf2013_official_residual_terms.py",
            "derived_outputs": "MF2013 basic and D1400/AVS30 station residuals",
            "access_date": TODAY,
            "redistribution_note": "Code and derived residual summaries can be released; cite the original model.",
        },
        {
            "asset": "J-SHIS hazard-curve API",
            "role": "Official hazard-curve reference",
            "version_or_date": "Y2024 AVR TTL_MTTL API",
            "source": "J-SHIS map API",
            "source_url_or_identifier": "https://www.j-shis.bosai.go.jp/map/api/pshm/Y2024/AVR/TTL_MTTL/T50/hzcv.json",
            "local_reference": "work/external_data/jshis_hazard_api_cache/",
            "derived_outputs": "jshis_official_hazard_curve_api_examples.csv",
            "access_date": TODAY,
            "redistribution_note": "Cache is a reproducibility aid; cite and refetch from J-SHIS for final release if required.",
        },
        {
            "asset": "J-SHIS response-spectrum map",
            "role": "Official SA spectrum ordinates for station correction check",
            "version_or_date": "Y2020 AVR TTL_MTTL T50 response-map ZIP",
            "source": "J-SHIS response-spectrum map download",
            "source_url_or_identifier": "https://www.j-shis.bosai.go.jp/",
            "local_reference": "work/external_data/jshis_respmap/P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip",
            "derived_outputs": "jshis_official_response_sa3_station_values.csv; jshis_official_hazard_response_summary.csv",
            "access_date": TODAY,
            "redistribution_note": "Release derived values and scripts; cite the J-SHIS product.",
        },
        {
            "asset": "Zhao et al. (2006) OpenQuake implementation",
            "role": "External GMPE sensitivity audit",
            "version_or_date": "OpenQuake hazardlib implementation used in existing project run",
            "source": "Zhao et al. (2006), Bulletin of the Seismological Society of America; OpenQuake hazardlib",
            "source_url_or_identifier": "https://doi.org/10.1785/0120050122",
            "local_reference": "outputs/jshis_zhao2006_external_gmpe_*.csv",
            "derived_outputs": "external-model residual metrics and station-site correlations",
            "access_date": TODAY,
            "redistribution_note": "Current environment cannot reimport OpenQuake because of geospatial-library requirements; existing derived outputs are retained with provenance.",
        },
        {
            "asset": "Station-residual model code and derived tables",
            "role": "Reproducible manuscript outputs",
            "version_or_date": "v0.8 manuscript package, generated 2026-06-14",
            "source": "This project repository",
            "source_url_or_identifier": "GitHub/Zenodo DOI to be filled after public release",
            "local_reference": "outputs/cee_submission_latex_v0_8_english_article/; outputs/nonergodic_station_model_repro_pack_v0_1.zip",
            "derived_outputs": "main figures, supplementary CSV tables, LaTeX source, reproducibility pack",
            "access_date": TODAY,
            "redistribution_note": "No DOI has been assigned yet; replace the placeholder after GitHub/Zenodo release.",
        },
    ]
    frame = pd.DataFrame(rows)
    frame.to_csv(output_dir / "jshis_official_data_provenance_2026-06-14.csv", index=False)
    lines = [
        "# Data and code release manifest",
        "",
        f"Generated on {TODAY}. This manifest separates original public products from derived project outputs. No DOI is claimed before deposition.",
        "",
        df_to_markdown(frame[["asset", "role", "version_or_date", "local_reference", "redistribution_note"]]),
        "",
        "Release action before submission: upload the code, derived CSV tables, LaTeX source, figures, and README to a public repository; archive the repository with Zenodo; replace the placeholder DOI in the manuscript and cover letter.",
        "",
    ]
    (output_dir / "cee_code_data_release_manifest_2026-06-14.md").write_text("\n".join(lines), encoding="utf-8")

    zenodo = {
        "metadata": {
            "title": "Station residual correction checks for Japanese strong-motion records and official hazard spectra",
            "upload_type": "dataset",
            "description": (
                "Derived station-level residual tables, official J-SHIS response-spectrum checks, "
                "and scripts supporting the manuscript. Original J-SHIS/NIED products should be "
                "obtained from the official provider."
            ),
            "creators": [
                {
                    "name": "Noda, Yojiro",
                    "affiliation": "To be confirmed before deposition",
                    "orcid": "To be confirmed",
                }
            ],
            "keywords": [
                "strong-motion records",
                "station residuals",
                "non-ergodic ground motion",
                "Japan",
                "J-SHIS",
                "hazard spectra",
            ],
            "license": "cc-by-4.0",
            "access_right": "open",
            "related_identifiers": [
                {
                    "identifier": "https://www.j-shis.bosai.go.jp/",
                    "relation": "isDerivedFrom",
                    "resource_type": "dataset",
                },
                {
                    "identifier": "https://doi.org/10.1785/0120120097",
                    "relation": "cites",
                    "resource_type": "publication-article",
                },
                {
                    "identifier": "https://doi.org/10.1785/0120050122",
                    "relation": "cites",
                    "resource_type": "publication-article",
                },
            ],
            "notes": "Replace affiliation, ORCID, repository URL, and DOI fields after author confirmation.",
        }
    }
    (output_dir / "cee_release_metadata_zenodo_draft.json").write_text(
        json.dumps(zenodo, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def run_existing_sensitivity_summary(output_dir: Path) -> None:
    lines = ["# Existing influence and external-model sensitivity summary", "", f"Generated on {TODAY}.", ""]
    tables = [
        ("Event-exclusion sensitivity", output_dir / "jshis_mf2013_event_exclusion_sensitivity.csv"),
        ("Station-exclusion sensitivity", output_dir / "jshis_mf2013_station_exclusion_sensitivity.csv"),
        ("Magnitude-bin sensitivity", output_dir / "jshis_mf2013_magnitude_bin_sensitivity.csv"),
        ("Distance-bin sensitivity", output_dir / "jshis_mf2013_distance_bin_sensitivity.csv"),
        ("Zhao et al. (2006) external GMPE overall metrics", output_dir / "jshis_zhao2006_external_gmpe_overall_metrics.csv"),
        ("Zhao et al. (2006) station-site correlations", output_dir / "jshis_zhao2006_external_gmpe_station_site_correlations.csv"),
    ]
    rows = []
    for label, path in tables:
        if not path.exists():
            rows.append({"check": label, "status": "missing", "rows": 0, "path": str(path.relative_to(PROJECT_ROOT))})
            continue
        data = pd.read_csv(path)
        rows.append({"check": label, "status": "available", "rows": int(data.shape[0]), "path": str(path.relative_to(PROJECT_ROOT))})
        preview = data.head(6)
        lines.extend([f"## {label}", "", df_to_markdown(preview), ""])
    pd.DataFrame(rows).to_csv(output_dir / "jshis_existing_sensitivity_inventory.csv", index=False)
    lines.extend(
        [
            "The current environment reports an OpenQuake geospatial-library import error, so no new external GMPE is generated here. The existing Zhao et al. (2006) derived outputs remain an external-model audit with explicit distance and VS30-proxy limitations.",
            "",
        ]
    )
    (output_dir / "jshis_existing_influence_and_external_model_sensitivity.md").write_text("\n".join(lines), encoding="utf-8")


def write_additional_experiment_report(output_dir: Path) -> None:
    region = pd.read_csv(output_dir / "jshis_station_region_leaveout_summary.csv")
    substitution = pd.read_csv(output_dir / "jshis_station_variable_substitution_summary.csv")
    positive = pd.read_csv(output_dir / "jshis_positive_sa3_station_correction_summary.csv")
    matching = pd.read_csv(output_dir / "jshis_station_matching_missingness_audit.csv")
    sa3_region = region[(region["target_label"].eq("SA(3.0s) RotD50")) & (region["model"].eq("full_site_location"))]
    sa3_official = substitution[
        (substitution["target_label"].eq("SA(3.0s) RotD50")) & (substitution["model"].isin(["mf2013_d1400_avs30", "official_four", "full_site_location"]))
    ]
    key_matching = matching[
        matching["check"].isin(["strong_motion_records", "unique_record_sites", "site_schema_rows", "record_site_schema_matched_rows"])
    ]
    lines = [
        "# Additional CEE submission experiments",
        "",
        f"Generated on {TODAY}. These tables cover the required, recommended, and feasible optional checks requested before CEE submission.",
        "",
        "## Regional extrapolation",
        "",
        df_to_markdown(sa3_region),
        "",
        "## Official-variable substitution for SA(3.0 s)",
        "",
        df_to_markdown(sa3_official),
        "",
        "## Positive SA(3.0 s) corrections",
        "",
        df_to_markdown(positive.head(12)),
        "",
        "## Matching audit",
        "",
        df_to_markdown(key_matching),
        "",
        "Full CSV outputs contain all periods and all feature sets.",
        "",
    ]
    (output_dir / "jshis_cee_additional_experiments_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--station-terms", type=Path, default=DEFAULT_STATION_TERMS)
    parser.add_argument("--site-schema", type=Path, default=DEFAULT_SITE_SCHEMA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seed", type=int, default=20260614)
    parser.add_argument("--min-records", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = load_station_dataset(args.station_terms, args.site_schema, min_records=args.min_records)
    frame = add_region(frame)
    run_region_leaveout(frame, args.output_dir, args.seed)
    run_variable_substitution(frame, args.output_dir, args.seed)
    run_collinearity(frame, args.output_dir)
    run_positive_correction_audit(frame, args.output_dir)
    run_matching_missingness_audit(args.output_dir, args.site_schema, args.station_terms)
    run_provenance_release_manifest(args.output_dir)
    run_existing_sensitivity_summary(args.output_dir)
    write_additional_experiment_report(args.output_dir)
    print(f"Additional CEE checks written to {args.output_dir}")


if __name__ == "__main__":
    main()
