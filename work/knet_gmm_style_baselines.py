#!/usr/bin/env python3
"""GMM-style baselines for K-NET waveform and response-spectrum targets.

This script implements an empirical ground-motion-model functional form using
magnitude, magnitude squared, hypocentral distance, log distance, depth, and
magnitude-distance interaction terms. It is not a replication of any published
GMPE/GMM. It is a stronger conventional baseline than the simple
magnitude-distance-depth ridge model used in earlier exploratory scripts.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_WAVEFORM_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_record_waveform_features.csv"
)
DEFAULT_SA_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_record_response_spectrum_features.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)

TARGETS = [
    ("log10_pga_horizontal_from_waveform_gal", "Horizontal PGA"),
    ("log10_arias_proxy_horizontal", "Horizontal Arias proxy"),
    ("log10_spectral_amp_1_2_hz_horizontal", "1-2 Hz Fourier amp"),
    ("log10_psa_0.3s_5pct_horizontal", "SA(0.3s)"),
    ("log10_psa_1s_5pct_horizontal", "SA(1.0s)"),
    ("log10_psa_3s_5pct_horizontal", "SA(3.0s)"),
]


def load_frame(waveform_path: Path, sa_path: Path) -> pd.DataFrame:
    wave = pd.read_csv(waveform_path, dtype={"_id": str, "EventName": str})
    sa_cols = [
        "_id",
        "log10_psa_0.3s_5pct_horizontal",
        "log10_psa_1s_5pct_horizontal",
        "log10_psa_3s_5pct_horizontal",
    ]
    sa = pd.read_csv(sa_path, usecols=sa_cols, dtype={"_id": str})
    df = wave.merge(sa, on="_id", how="left")
    df = df[
        df["Magnitude"].notna()
        & df["Distance"].notna()
        & df["Depth_km"].notna()
        & df["StationCode"].notna()
        & df["EventName"].notna()
        & df["Latitude"].notna()
        & df["Longitude"].notna()
        & df["StationLat"].notna()
        & df["StationLong"].notna()
    ].copy()

    m = df["Magnitude"].astype(float)
    r_epi = df["Distance"].astype(float)
    depth = df["Depth_km"].astype(float)
    r_hyp = np.sqrt(r_epi**2 + depth**2)
    log_r = np.log10(r_hyp + 1.0)
    df["m_centered"] = m - 6.0
    df["m_centered_sq"] = df["m_centered"] ** 2
    df["hyp_distance_km"] = r_hyp
    df["log10_hyp_distance_km"] = log_r
    df["linear_distance_100km"] = r_hyp / 100.0
    df["depth_50km"] = depth / 50.0
    df["m_log_distance_interaction"] = df["m_centered"] * log_r
    df["shallow_depth_flag"] = (depth <= 30.0).astype(float)
    return df


def pipeline(numeric: list[str], categorical: list[str] | None = None, alpha: float = 1.0) -> Pipeline:
    categorical = categorical or []
    transformers = [
        (
            "num",
            Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]),
            numeric,
        )
    ]
    if categorical:
        transformers.append(
            (
                "cat",
                OneHotEncoder(handle_unknown="ignore", min_frequency=2),
                categorical,
            )
        )
    return Pipeline(
        [
            ("features", ColumnTransformer(transformers=transformers, remainder="drop")),
            ("model", Ridge(alpha=alpha)),
        ]
    )


def hgb_pipeline(numeric: list[str]) -> Pipeline:
    return Pipeline(
        [
            (
                "features",
                ColumnTransformer(
                    [("num", Pipeline([("impute", SimpleImputer(strategy="median"))]), numeric)],
                    remainder="drop",
                ),
            ),
            (
                "model",
                HistGradientBoostingRegressor(
                    max_iter=300,
                    learning_rate=0.05,
                    max_leaf_nodes=31,
                    l2_regularization=0.05,
                    random_state=20260605,
                ),
            ),
        ]
    )


def build_models() -> dict[str, Pipeline]:
    gmm = [
        "m_centered",
        "m_centered_sq",
        "log10_hyp_distance_km",
        "linear_distance_100km",
        "depth_50km",
        "m_log_distance_interaction",
        "shallow_depth_flag",
    ]
    geo = gmm + [
        "Latitude",
        "Longitude",
        "StationLat",
        "StationLong",
        "StationHeight_m",
    ]
    return {
        "gmm_style": pipeline(gmm, alpha=0.1),
        "gmm_style_site_terms": pipeline(geo, ["StationCode"], alpha=5.0),
        "hgb_geo_path": hgb_pipeline(geo),
    }


def split_indices(df: pd.DataFrame, split_col: str) -> tuple[np.ndarray, np.ndarray]:
    train = df[split_col].isin(["train", "dev"]).to_numpy()
    test = (df[split_col] == "test").to_numpy()
    return np.where(train)[0], np.where(test)[0]


def evaluate(df: pd.DataFrame, split: str, target: str, label: str, model_name: str, model: Pipeline) -> tuple[dict, pd.DataFrame]:
    work = df[df[target].notna()].copy()
    train_idx, test_idx = split_indices(work, split)
    y = work[target].to_numpy()
    model.fit(work.iloc[train_idx], y[train_idx])
    pred = model.predict(work.iloc[test_idx])
    metrics = {
        "split": split,
        "target": target,
        "target_label": label,
        "model": model_name,
        "train_rows": int(len(train_idx)),
        "test_rows": int(len(test_idx)),
        "mae": float(mean_absolute_error(y[test_idx], pred)),
        "r2": float(r2_score(y[test_idx], pred)),
    }
    pred_df = work.iloc[test_idx][["EventName", "StationCode"]].copy()
    pred_df["split"] = split
    pred_df["target"] = target
    pred_df["target_label"] = label
    pred_df["model"] = model_name
    pred_df["absolute_error"] = np.abs(y[test_idx] - pred)
    return metrics, pred_df


def add_reductions(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split, target), sub in results.groupby(["split", "target"]):
        base = sub[sub["model"] == "gmm_style"].iloc[0]
        for _, row in sub.iterrows():
            item = row.to_dict()
            item["mae_reduction_vs_gmm_style"] = base["mae"] - row["mae"]
            item["percent_mae_reduction_vs_gmm_style"] = (
                (base["mae"] - row["mae"]) / base["mae"] * 100.0
            )
            rows.append(item)
    return pd.DataFrame(rows)


def bootstrap_ci(predictions: pd.DataFrame, n_bootstrap: int = 400) -> pd.DataFrame:
    rng = np.random.default_rng(20260605)
    rows = []
    for (split, target, model), sub in predictions.groupby(["split", "target", "model"]):
        group_col = "StationCode" if split == "station_split" else "EventName"
        groups = sub[group_col].astype(str).unique()
        errors_by_group = {
            group: sub.loc[sub[group_col].astype(str) == group, "absolute_error"].to_numpy()
            for group in groups
        }
        vals = []
        for _ in range(n_bootstrap):
            sampled = rng.choice(groups, size=len(groups), replace=True)
            errs = np.concatenate([errors_by_group[group] for group in sampled])
            vals.append(float(np.mean(errs)))
        vals = np.asarray(vals)
        rows.append(
            {
                "split": split,
                "target": target,
                "model": model,
                "group_col": group_col,
                "n_groups": int(len(groups)),
                "mae_ci_low": float(np.quantile(vals, 0.025)),
                "mae_ci_high": float(np.quantile(vals, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def plot_heatmaps(reductions: pd.DataFrame, output_dir: Path) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for split in ["event_split", "time_split", "station_split"]:
        sub = reductions[(reductions["split"] == split) & (reductions["model"] != "gmm_style")]
        pivot = sub.pivot(
            index="target_label",
            columns="model",
            values="percent_mae_reduction_vs_gmm_style",
        )
        plt.figure(figsize=(7.4, 4.8))
        im = plt.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
        plt.yticks(np.arange(len(pivot.index)), pivot.index)
        plt.xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=20, ha="right")
        plt.colorbar(im, label="MAE reduction vs GMM-style baseline (%)")
        plt.title(f"GMM-style baseline gains: {split}")
        plt.tight_layout()
        plt.savefig(figures / f"knet_gmm_style_reduction_heatmap_{split}.png", dpi=220)
        plt.close()


def write_report(reductions: pd.DataFrame, ci: pd.DataFrame, output_dir: Path) -> None:
    merged = reductions.merge(ci, on=["split", "target", "model"], how="left")
    lines = ["# K-NET GMM-Style Baseline Comparison\n\n"]
    lines.append(
        "The `gmm_style` baseline is an empirical functional form with magnitude, magnitude squared, log hypocentral distance, linear distance, depth, and magnitude-distance interaction terms. It is not a published GMPE replication.\n\n"
    )
    for split in ["event_split", "time_split", "station_split"]:
        lines.append(f"## {split}\n")
        sub = merged[merged["split"] == split]
        for _, label in TARGETS:
            rows = sub[sub["target_label"] == label].sort_values("mae")
            best = rows.iloc[0]
            base = rows[rows["model"] == "gmm_style"].iloc[0]
            lines.append(
                f"- {label}: best={best['model']} MAE={best['mae']:.4f}, "
                f"95% CI={best['mae_ci_low']:.4f}-{best['mae_ci_high']:.4f}, "
                f"R2={best['r2']:.4f}, gain vs GMM-style={best['percent_mae_reduction_vs_gmm_style']:.1f}% "
                f"(GMM-style MAE={base['mae']:.4f})\n"
            )
        lines.append("\n")
    lines.append("## Guardrails\n")
    lines.append("- This is a fitted empirical GMM-style baseline, not an external published GMPE.\n")
    lines.append("- The next step is to implement or cite a published Japan/NIED-compatible GMPE for comparison.\n")
    lines.append("- Gains over GMM-style are more defensible than gains over the earlier simple ridge physical baseline.\n")
    (output_dir / "knet_gmm_style_baseline_comparison.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--waveform-input", type=Path, default=DEFAULT_WAVEFORM_INPUT)
    parser.add_argument("--sa-input", type=Path, default=DEFAULT_SA_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_frame(args.waveform_input, args.sa_input)
    models = build_models()
    rows = []
    pred_rows = []
    for split in ["event_split", "time_split", "station_split"]:
        for target, label in TARGETS:
            for model_name, model in models.items():
                metrics, preds = evaluate(df, split, target, label, model_name, model)
                rows.append(metrics)
                pred_rows.append(preds)
    results = pd.DataFrame(rows)
    predictions = pd.concat(pred_rows, ignore_index=True)
    reductions = add_reductions(results)
    ci = bootstrap_ci(predictions)
    results.to_csv(args.output_dir / "knet_gmm_style_baseline_metrics.csv", index=False)
    reductions.to_csv(args.output_dir / "knet_gmm_style_baseline_reductions.csv", index=False)
    ci.to_csv(args.output_dir / "knet_gmm_style_baseline_bootstrap_ci.csv", index=False)
    plot_heatmaps(reductions, args.output_dir)
    write_report(reductions, ci, args.output_dir)


if __name__ == "__main__":
    main()
