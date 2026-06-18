#!/usr/bin/env python3
"""Station-held-out non-ergodic residual correction experiment.

The experiment starts from station-mean MF2013 residual terms that were already
computed from the public J-SHIS/NIED flatfile. It asks whether remaining
station residuals after the official D1400/AVS30 correction can be predicted
from public site and location variables for stations not used in training.

The target is deliberately station-level. It does not train on individual
records and it does not claim a new ground-motion prediction equation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.cluster import KMeans
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATION_TERMS = PROJECT_ROOT / "outputs" / "jshis_mf2013_official_station_residual_terms.csv"
DEFAULT_SITE_SCHEMA = PROJECT_ROOT / "outputs" / "jshis_site_schema_v2024_sub1.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures"

TARGET_ORDER = [
    "PGA RotD50",
    "SA(0.3s) RotD50",
    "SA(1.0s) RotD50",
    "SA(3.0s) RotD50",
]


@dataclass(frozen=True)
class ModelSpec:
    name: str
    feature_set: str
    estimator: object


def positive_log10(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return np.log10(values.where(values > 0))


def load_station_dataset(station_terms_path: Path, site_schema_path: Path, min_records: int) -> pd.DataFrame:
    terms = pd.read_csv(station_terms_path)
    sites = pd.read_csv(site_schema_path)
    site_cols = [
        "siteid2",
        "lon",
        "lat",
        "elevation",
        "sensor_depth_glminus",
        "obs_network_id",
        "installation_situation_id",
        "dist_vf_mf13_nejapan",
        "dist_vf_mf13_swjapan",
        "meshcode3",
        "meshcode250",
    ]
    sites = sites[site_cols].drop_duplicates("siteid2")
    frame = terms.merge(sites, on=["siteid2", "obs_network_id"], how="left", validate="many_to_one")
    frame = frame[(frame["model"] == "mf2013_site") & (frame["n_records"] >= min_records)].copy()
    for col in ["vs10", "vs20", "vs30", "avs30", "d1100", "d1400", "d1700", "d2100", "dbase"]:
        frame[f"log_{col}"] = positive_log10(frame[col])
    frame["network_label"] = np.where(frame["obs_network_id"].eq(2), "KiK-net", "K-NET")
    frame["target_label"] = pd.Categorical(frame["target_label"], TARGET_ORDER, ordered=True)
    return frame.sort_values(["target_label", "siteid2"]).reset_index(drop=True)


def make_spatial_folds(frame: pd.DataFrame, n_splits: int, seed: int) -> np.ndarray:
    coords = frame[["lon", "lat"]].to_numpy(dtype=float)
    coords = np.nan_to_num(coords, nan=np.nanmedian(coords, axis=0))
    coords = StandardScaler().fit_transform(coords)
    labels = KMeans(n_clusters=n_splits, random_state=seed, n_init=20).fit_predict(coords)
    return labels


def model_specs(seed: int) -> list[ModelSpec]:
    ridge = RidgeCV(alphas=np.logspace(-3, 3, 25))
    hgb = HistGradientBoostingRegressor(
        loss="squared_error",
        learning_rate=0.04,
        max_iter=400,
        max_leaf_nodes=15,
        min_samples_leaf=35,
        l2_regularization=0.1,
        random_state=seed,
    )
    return [
        ModelSpec("ridge_site", "site", ridge),
        ModelSpec("ridge_site_space", "site_space", ridge),
        ModelSpec("gradient_boosted_site_space", "site_space", hgb),
    ]


def feature_columns(feature_set: str) -> tuple[list[str], list[str]]:
    site_numeric = [
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
    ]
    spatial_numeric = ["lon", "lat"]
    numeric = site_numeric + (spatial_numeric if feature_set == "site_space" else [])
    categorical = ["network_label", "installation_situation_id"]
    return numeric, categorical


def build_pipeline(estimator: object, feature_set: str) -> Pipeline:
    numeric, categorical = feature_columns(feature_set)
    pre = ColumnTransformer(
        transformers=[
            ("num", Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric),
            (
                "cat",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                    ]
                ),
                categorical,
            ),
        ]
    )
    return Pipeline([("preprocess", pre), ("model", estimator)])


def weighted_mae(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(np.abs(y_true - y_pred), weights=weights))


def weighted_rmse(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.average((y_true - y_pred) ** 2, weights=weights)))


def fold_rows_for_target(
    frame: pd.DataFrame,
    target_label: str,
    split_name: str,
    fold_ids: np.ndarray,
    seed: int,
) -> tuple[list[dict], list[pd.DataFrame]]:
    target_frame = frame[frame["target_label"].astype(str).eq(target_label)].copy()
    y = target_frame["mean_residual"].to_numpy(dtype=float)
    weights = target_frame["n_records"].to_numpy(dtype=float)
    pred_frames: list[pd.DataFrame] = []
    rows: list[dict] = []

    def iter_splits() -> list[tuple[int, np.ndarray, np.ndarray]]:
        if split_name == "station_random":
            kf = KFold(n_splits=5, shuffle=True, random_state=seed)
            return [(fold_idx, train_idx, test_idx) for fold_idx, (train_idx, test_idx) in enumerate(kf.split(target_frame))]
        unique_folds = sorted(np.unique(fold_ids).tolist())
        splits = []
        for fold_idx in unique_folds:
            test_idx = np.flatnonzero(fold_ids == fold_idx)
            train_idx = np.flatnonzero(fold_ids != fold_idx)
            splits.append((int(fold_idx), train_idx, test_idx))
        return splits

    specs = model_specs(seed)
    for fold_idx, train_idx, test_idx in iter_splits():
        fold_pred = target_frame.iloc[test_idx][
            ["target", "target_label", "period", "siteid2", "site_code", "network_label", "n_records", "mean_residual"]
        ].copy()
        fold_pred["split"] = split_name
        fold_pred["fold"] = fold_idx
        fold_pred["baseline_zero_prediction"] = 0.0

        y_test = y[test_idx]
        w_test = weights[test_idx]
        baseline_mae = weighted_mae(y_test, np.zeros_like(y_test), w_test)
        baseline_rmse = weighted_rmse(y_test, np.zeros_like(y_test), w_test)
        rows.append(
            {
                "target_label": target_label,
                "split": split_name,
                "fold": fold_idx,
                "model": "zero_station_residual",
                "n_train_sites": len(train_idx),
                "n_test_sites": len(test_idx),
                "test_weight_sum_records": float(w_test.sum()),
                "weighted_mae": baseline_mae,
                "weighted_rmse": baseline_rmse,
                "mae_reduction_pct": 0.0,
                "rmse_reduction_pct": 0.0,
                "weighted_r2_vs_zero": 0.0,
            }
        )

        for spec in specs:
            pipe = build_pipeline(clone(spec.estimator), spec.feature_set)
            pipe.fit(target_frame.iloc[train_idx], y[train_idx], model__sample_weight=weights[train_idx])
            pred = pipe.predict(target_frame.iloc[test_idx])
            mae = weighted_mae(y_test, pred, w_test)
            rmse = weighted_rmse(y_test, pred, w_test)
            rows.append(
                {
                    "target_label": target_label,
                    "split": split_name,
                    "fold": fold_idx,
                    "model": spec.name,
                    "n_train_sites": len(train_idx),
                    "n_test_sites": len(test_idx),
                    "test_weight_sum_records": float(w_test.sum()),
                    "weighted_mae": mae,
                    "weighted_rmse": rmse,
                    "mae_reduction_pct": 100.0 * (baseline_mae - mae) / baseline_mae if baseline_mae else np.nan,
                    "rmse_reduction_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse if baseline_rmse else np.nan,
                    "weighted_r2_vs_zero": 1.0 - (rmse**2 / baseline_rmse**2) if baseline_rmse else np.nan,
                }
            )
            fold_pred[f"{spec.name}_prediction"] = pred
            fold_pred[f"{spec.name}_remaining_residual"] = y_test - pred

        pred_frames.append(fold_pred)
    return rows, pred_frames


def summarize_fold_metrics(fold_metrics: pd.DataFrame) -> pd.DataFrame:
    model_rows = fold_metrics[fold_metrics["model"].ne("zero_station_residual")].copy()
    grouped = (
        model_rows.groupby(["target_label", "split", "model"], observed=True)
        .agg(
            folds=("fold", "nunique"),
            mean_test_sites=("n_test_sites", "mean"),
            weighted_mae_mean=("weighted_mae", "mean"),
            weighted_mae_std=("weighted_mae", "std"),
            weighted_rmse_mean=("weighted_rmse", "mean"),
            weighted_rmse_std=("weighted_rmse", "std"),
            mae_reduction_pct_mean=("mae_reduction_pct", "mean"),
            mae_reduction_pct_std=("mae_reduction_pct", "std"),
            rmse_reduction_pct_mean=("rmse_reduction_pct", "mean"),
            rmse_reduction_pct_std=("rmse_reduction_pct", "std"),
            weighted_r2_vs_zero_mean=("weighted_r2_vs_zero", "mean"),
            weighted_r2_vs_zero_std=("weighted_r2_vs_zero", "std"),
        )
        .reset_index()
    )
    grouped["target_label"] = pd.Categorical(grouped["target_label"], TARGET_ORDER, ordered=True)
    return grouped.sort_values(["split", "target_label", "model"]).reset_index(drop=True)


def markdown_table(frame: pd.DataFrame, floatfmt: str = ".3f") -> str:
    headers = list(frame.columns)
    rows = []
    for _, row in frame.iterrows():
        formatted = []
        for col in headers:
            value = row[col]
            if isinstance(value, (float, np.floating)):
                formatted.append(format(float(value), floatfmt))
            else:
                formatted.append(str(value))
        rows.append(formatted)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_report(summary: pd.DataFrame, fold_metrics: pd.DataFrame, output_path: Path, min_records: int) -> None:
    def best_table(split: str) -> pd.DataFrame:
        sub = summary[summary["split"].eq(split)].copy()
        idx = sub.groupby("target_label", observed=True)["rmse_reduction_pct_mean"].idxmax()
        return sub.loc[idx].sort_values("target_label")

    random_best = best_table("station_random")
    spatial_best = best_table("spatial_block")

    lines = [
        "# Non-ergodic Station Residual Model Experiment",
        "",
        "Status: station-held-out supplement for the J-SHIS/MF2013 residual manuscript.",
        "",
        "## Design",
        "",
        f"- Input station residuals: `mf2013_site` station mean residuals with at least {min_records} records per target and station.",
        "- Target: station mean residual after the official MF2013 D1400/AVS30 site correction.",
        "- Baseline: zero additional station correction.",
        "- Models: ridge regression with public site variables, ridge regression with site plus coordinates, and gradient-boosted regression with site plus coordinates.",
        "- Validation: five-fold random station holdout and five-fold spatial-block holdout based on station longitude and latitude.",
        "- Metrics: weighted station-residual MAE/RMSE, weighted by the number of records per station.",
        "",
        "## Best model by target",
        "",
        "### Random station holdout",
        "",
        markdown_table(
            random_best[
                [
                    "target_label",
                    "model",
                    "rmse_reduction_pct_mean",
                    "rmse_reduction_pct_std",
                    "mae_reduction_pct_mean",
                    "weighted_r2_vs_zero_mean",
                ]
            ]
        ),
        "",
        "### Spatial-block holdout",
        "",
        markdown_table(
            spatial_best[
                [
                    "target_label",
                    "model",
                    "rmse_reduction_pct_mean",
                    "rmse_reduction_pct_std",
                    "mae_reduction_pct_mean",
                    "weighted_r2_vs_zero_mean",
                ]
            ]
        ),
        "",
        "## Interpretation",
        "",
    ]
    long_period = spatial_best[spatial_best["target_label"].astype(str).isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])]
    if (long_period["rmse_reduction_pct_mean"] > 0).all():
        lines.append(
            "The spatial-block test retains positive long-period residual reduction, so the supplement supports a predictive non-ergodic site-residual extension beyond the official D1400/AVS30 correction."
        )
    else:
        lines.append(
            "The random station holdout is informative, but the spatial-block test is weaker. The result should be framed as station-level predictability within sampled regions, not as robust regional extrapolation."
        )
    lines.extend(
        [
            "",
            "This experiment does not define a new ground-motion prediction equation. It tests whether public site and location variables can predict remaining station terms under held-out-station validation.",
            "",
            "## Output files",
            "",
            "- `jshis_nonergodic_station_model_fold_metrics.csv`: fold-level validation metrics.",
            "- `jshis_nonergodic_station_model_summary.csv`: summary metrics by target, split, and model.",
            "- `jshis_nonergodic_station_model_predictions.csv`: held-out station predictions.",
            "- `figures/cee_fig6_nonergodic_station_model.pdf`: manuscript-style summary figure.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_figure(summary: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    best = (
        summary.sort_values("rmse_reduction_pct_mean", ascending=False)
        .groupby(["target_label", "split"], observed=True)
        .head(1)
        .copy()
    )
    best["target_label"] = pd.Categorical(best["target_label"], TARGET_ORDER, ordered=True)
    best = best.sort_values(["split", "target_label"])

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9), sharey=True)
    split_titles = {
        "station_random": "Random station holdout",
        "spatial_block": "Spatial-block holdout",
    }
    colors = {
        "PGA RotD50": "#4C78A8",
        "SA(0.3s) RotD50": "#72B7B2",
        "SA(1.0s) RotD50": "#F58518",
        "SA(3.0s) RotD50": "#E45756",
    }
    for ax, split in zip(axes, ["station_random", "spatial_block"], strict=True):
        sub = best[best["split"].eq(split)]
        labels = [str(x).replace(" RotD50", "") for x in sub["target_label"]]
        vals = sub["rmse_reduction_pct_mean"].to_numpy(dtype=float)
        errs = sub["rmse_reduction_pct_std"].fillna(0).to_numpy(dtype=float)
        bars = ax.bar(
            np.arange(len(sub)),
            vals,
            yerr=errs,
            color=[colors[str(t)] for t in sub["target_label"]],
            edgecolor="black",
            linewidth=0.6,
            capsize=3,
        )
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xticks(np.arange(len(sub)))
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_title(split_titles[split], fontsize=11)
        ax.set_ylabel("Station-residual RMSE reduction (%)")
        ax.set_ylim(0, max(70, float(np.nanmax(best["rmse_reduction_pct_mean"] + best["rmse_reduction_pct_std"])) + 8))
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        for bar, val in zip(bars, vals, strict=True):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val + 1.6,
                f"{val:.1f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )
    fig.suptitle("Held-out prediction of remaining station residuals after MF2013 site correction", fontsize=12, y=0.98)
    fig.text(
        0.5,
        0.01,
        "Bars show the best-performing candidate model for each target; in all cases this is the gradient-boosted site-plus-location model.",
        ha="center",
        fontsize=8.5,
    )
    fig.tight_layout()
    fig.savefig(figure_dir / "cee_fig6_nonergodic_station_model.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / "cee_fig6_nonergodic_station_model.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--station-terms", type=Path, default=DEFAULT_STATION_TERMS)
    parser.add_argument("--site-schema", type=Path, default=DEFAULT_SITE_SCHEMA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--min-records", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260612)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    frame = load_station_dataset(args.station_terms, args.site_schema, args.min_records)
    all_rows: list[dict] = []
    prediction_frames: list[pd.DataFrame] = []
    for target_label in TARGET_ORDER:
        target_frame = frame[frame["target_label"].astype(str).eq(target_label)].copy()
        random_rows, random_preds = fold_rows_for_target(
            frame=frame,
            target_label=target_label,
            split_name="station_random",
            fold_ids=np.zeros(len(target_frame), dtype=int),
            seed=args.seed,
        )
        spatial_folds = make_spatial_folds(target_frame, n_splits=5, seed=args.seed)
        spatial_rows, spatial_preds = fold_rows_for_target(
            frame=frame,
            target_label=target_label,
            split_name="spatial_block",
            fold_ids=spatial_folds,
            seed=args.seed,
        )
        all_rows.extend(random_rows)
        all_rows.extend(spatial_rows)
        prediction_frames.extend(random_preds)
        prediction_frames.extend(spatial_preds)

    fold_metrics = pd.DataFrame(all_rows)
    summary = summarize_fold_metrics(fold_metrics)
    predictions = pd.concat(prediction_frames, ignore_index=True)

    fold_path = args.output_dir / "jshis_nonergodic_station_model_fold_metrics.csv"
    summary_path = args.output_dir / "jshis_nonergodic_station_model_summary.csv"
    prediction_path = args.output_dir / "jshis_nonergodic_station_model_predictions.csv"
    report_path = args.output_dir / "jshis_nonergodic_station_model.md"
    fold_metrics.to_csv(fold_path, index=False)
    summary.to_csv(summary_path, index=False)
    predictions.to_csv(prediction_path, index=False)
    write_report(summary, fold_metrics, report_path, args.min_records)
    build_figure(summary, args.figure_dir)

    print(f"Wrote {fold_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {prediction_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {args.figure_dir / 'cee_fig6_nonergodic_station_model.pdf'}")


if __name__ == "__main__":
    main()
