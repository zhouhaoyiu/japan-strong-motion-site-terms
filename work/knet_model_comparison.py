#!/usr/bin/env python3
"""Leakage-safe baseline comparisons for K-NET log10(PGA).

The goal is not to claim a final model. The goal is to identify whether source,
path, and site information adds measurable value under event/time/station
holdouts. These results define the next paper-grade experiments.
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


DEFAULT_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_header_index.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def load_frame(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["OriginTime", "RecordTime", "pTime", "sTime"])
    df = df[
        df["PGA_gal"].notna()
        & df["Magnitude"].notna()
        & df["Distance"].notna()
        & df["Depth_km"].notna()
        & df["StationCode"].notna()
        & df["EventName"].notna()
        & (df["PGA_gal"] > 0)
        & (df["Distance"] >= 0)
    ].copy()
    df["log10_pga_gal"] = np.log10(df["PGA_gal"].astype(float))
    df["log10_distance_km"] = np.log10(df["Distance"].astype(float) + 1.0)
    df["hyp_distance_km"] = np.sqrt(
        df["Distance"].astype(float) ** 2 + df["Depth_km"].astype(float) ** 2
    )
    df["log10_hyp_distance_km"] = np.log10(df["hyp_distance_km"] + 1.0)
    return df


def ridge_pipeline(numeric: list[str], categorical: list[str] | None = None) -> Pipeline:
    categorical = categorical or []
    transformers = [
        (
            "num",
            Pipeline(
                [
                    ("impute", SimpleImputer(strategy="median")),
                    ("scale", StandardScaler()),
                ]
            ),
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
    pre = ColumnTransformer(transformers=transformers, remainder="drop")
    return Pipeline([("features", pre), ("model", Ridge(alpha=5.0))])


def hist_pipeline(numeric: list[str]) -> Pipeline:
    pre = ColumnTransformer(
        transformers=[
            (
                "num",
                Pipeline([("impute", SimpleImputer(strategy="median"))]),
                numeric,
            )
        ],
        remainder="drop",
    )
    return Pipeline(
        [
            ("features", pre),
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


def split_indices(df: pd.DataFrame, split_col: str) -> tuple[np.ndarray, np.ndarray]:
    train_mask = df[split_col].isin(["train", "dev"]).to_numpy()
    test_mask = (df[split_col] == "test").to_numpy()
    return np.where(train_mask)[0], np.where(test_mask)[0]


def evaluate_model(
    df: pd.DataFrame, split_col: str, model_name: str, model: Pipeline
) -> tuple[dict, pd.DataFrame]:
    train_idx, test_idx = split_indices(df, split_col)
    y = df["log10_pga_gal"].to_numpy()
    model.fit(df.iloc[train_idx], y[train_idx])
    pred = model.predict(df.iloc[test_idx])
    metrics = {
        "split": split_col,
        "model": model_name,
        "train_rows": int(len(train_idx)),
        "test_rows": int(len(test_idx)),
        "mae_log10_pga": float(mean_absolute_error(y[test_idx], pred)),
        "r2_log10_pga": float(r2_score(y[test_idx], pred)),
    }
    pred_df = df.iloc[test_idx][["EventName", "StationCode", "time_split", "event_split", "station_split"]].copy()
    pred_df["split"] = split_col
    pred_df["model"] = model_name
    pred_df["y_true_log10_pga"] = y[test_idx]
    pred_df["y_pred_log10_pga"] = pred
    pred_df["absolute_error_log10_pga"] = np.abs(pred_df["y_true_log10_pga"] - pred_df["y_pred_log10_pga"])
    return metrics, pred_df


def bootstrap_mae_ci(predictions: pd.DataFrame, n_bootstrap: int = 500) -> pd.DataFrame:
    rng = np.random.default_rng(20260605)
    rows = []
    for (split, model), sub in predictions.groupby(["split", "model"]):
        group_col = "StationCode" if split == "station_split" else "EventName"
        groups = sub[group_col].astype(str).unique()
        group_to_errors = {
            group: sub.loc[sub[group_col].astype(str) == group, "absolute_error_log10_pga"].to_numpy()
            for group in groups
        }
        vals = []
        for _ in range(n_bootstrap):
            sampled = rng.choice(groups, size=len(groups), replace=True)
            errs = np.concatenate([group_to_errors[group] for group in sampled])
            vals.append(float(np.mean(errs)))
        vals = np.array(vals)
        rows.append(
            {
                "split": split,
                "model": model,
                "group_col": group_col,
                "n_groups": int(len(groups)),
                "mae_bootstrap_mean": float(vals.mean()),
                "mae_ci_low": float(np.quantile(vals, 0.025)),
                "mae_ci_high": float(np.quantile(vals, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def baseline_reduction_table(results: pd.DataFrame, ci: pd.DataFrame) -> pd.DataFrame:
    merged = results.merge(ci, on=["split", "model"], how="left")
    rows = []
    for split, sub in merged.groupby("split"):
        baseline = sub[sub["model"] == "ridge_physical"].iloc[0]
        for _, row in sub.iterrows():
            rows.append(
                {
                    "split": split,
                    "model": row["model"],
                    "mae_log10_pga": row["mae_log10_pga"],
                    "mae_ci_low": row["mae_ci_low"],
                    "mae_ci_high": row["mae_ci_high"],
                    "absolute_mae_reduction_vs_physical": baseline["mae_log10_pga"]
                    - row["mae_log10_pga"],
                    "percent_mae_reduction_vs_physical": (
                        (baseline["mae_log10_pga"] - row["mae_log10_pga"])
                        / baseline["mae_log10_pga"]
                        * 100.0
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_models() -> dict[str, Pipeline]:
    physical = ["Magnitude", "log10_distance_km", "Depth_km"]
    physical_hyp = ["Magnitude", "log10_distance_km", "log10_hyp_distance_km", "Depth_km"]
    geo = physical_hyp + ["Latitude", "Longitude", "StationLat", "StationLong", "StationHeight_m"]
    return {
        "ridge_physical": ridge_pipeline(physical),
        "ridge_physical_hyp": ridge_pipeline(physical_hyp),
        "ridge_geo_path": ridge_pipeline(geo),
        "ridge_site_terms": ridge_pipeline(geo, ["StationCode"]),
        "hgb_geo_path": hist_pipeline(geo),
    }


def plot_results(results: pd.DataFrame, output_dir: Path) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    splits = ["event_split", "time_split", "station_split"]
    models = list(results["model"].drop_duplicates())
    x = np.arange(len(models))
    width = 0.24

    plt.figure(figsize=(10, 4.8))
    for offset, split in zip([-width, 0.0, width], splits):
        vals = [
            results[(results["split"] == split) & (results["model"] == model)][
                "mae_log10_pga"
            ].iloc[0]
            for model in models
        ]
        plt.bar(x + offset, vals, width=width, label=split)
    plt.xticks(x, models, rotation=25, ha="right")
    plt.ylabel("MAE in log10(PGA)")
    plt.title("Leakage-safe K-NET baseline comparison")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figures / "knet_model_comparison_mae.png", dpi=220)
    plt.close()


def write_summary(results: pd.DataFrame, reductions: pd.DataFrame, output_dir: Path) -> None:
    lines = ["# K-NET Leakage-Safe Model Comparison\n\n"]
    lines.append("Lower MAE in log10(PGA) is better. Test rows are disjoint by the named split.\n\n")
    for split in ["event_split", "time_split", "station_split"]:
        sub = results[results["split"] == split].sort_values("mae_log10_pga")
        lines.append(f"## {split}\n")
        for _, row in sub.iterrows():
            red = reductions[(reductions["split"] == split) & (reductions["model"] == row["model"])].iloc[0]
            lines.append(
                f"- {row['model']}: MAE={row['mae_log10_pga']:.4f}, "
                f"95% CI={red['mae_ci_low']:.4f}-{red['mae_ci_high']:.4f}, "
                f"R2={row['r2_log10_pga']:.4f}, "
                f"MAE reduction vs physical={red['percent_mae_reduction_vs_physical']:.1f}%, "
                f"test rows={int(row['test_rows']):,}\n"
            )
        lines.append("\n")
    lines.append("## Interpretation guardrails\n")
    lines.append(
        "- event_split tests generalization to unseen earthquakes; this is the primary predictive benchmark.\n"
    )
    lines.append(
        "- time_split tests historical-to-future transfer; this is closest to a publication-grade robustness check.\n"
    )
    lines.append(
        "- station_split tests unseen stations; station one-hot terms should not be over-interpreted under this split.\n"
    )
    lines.append(
        "- Better site-term performance under event/time splits is evidence for repeatable station effects, not yet proof of physical site amplification.\n"
    )
    (output_dir / "knet_model_comparison.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_frame(args.input)
    models = build_models()
    rows = []
    pred_rows = []
    for split_col in ["event_split", "time_split", "station_split"]:
        for model_name, model in models.items():
            metrics, pred_df = evaluate_model(df, split_col, model_name, model)
            rows.append(metrics)
            pred_rows.append(pred_df)
    results = pd.DataFrame(rows)
    results.to_csv(args.output_dir / "knet_model_comparison.csv", index=False)
    predictions = pd.concat(pred_rows, ignore_index=True)
    predictions.to_csv(args.output_dir / "knet_holdout_predictions.csv", index=False)
    ci = bootstrap_mae_ci(predictions)
    ci.to_csv(args.output_dir / "knet_model_comparison_bootstrap_ci.csv", index=False)
    reductions = baseline_reduction_table(results, ci)
    reductions.to_csv(args.output_dir / "knet_model_comparison_reductions.csv", index=False)
    plot_results(results, args.output_dir)
    write_summary(results, reductions, args.output_dir)


if __name__ == "__main__":
    main()
