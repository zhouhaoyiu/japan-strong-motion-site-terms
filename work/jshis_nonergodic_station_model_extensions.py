#!/usr/bin/env python3
"""Extensions for the station-held-out non-ergodic residual experiment.

This script adds targeted checks requested after the first held-out station
model:

1. feature-set ablation;
2. grouped permutation importance;
3. strict K-NET/KiK-net transfer validation;
4. engineering-scale residual correction examples;
5. a compact reproducibility pack for the new station-model outputs.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold

from jshis_nonergodic_station_residual_model import (
    DEFAULT_FIGURE_DIR,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_SITE_SCHEMA,
    DEFAULT_STATION_TERMS,
    TARGET_ORDER,
    build_pipeline,
    load_station_dataset,
    make_spatial_folds,
    weighted_mae,
    weighted_rmse,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = Path(__file__).resolve()
BASE_SCRIPT = PROJECT_ROOT / "work" / "jshis_nonergodic_station_residual_model.py"
HAZARD_SCRIPT = PROJECT_ROOT / "work" / "jshis_station_hazard_impact.py"
OFFICIAL_HAZARD_SCRIPT = PROJECT_ROOT / "work" / "jshis_official_hazard_response_check.py"
OFFICIAL_UNCERTAINTY_SCRIPT = PROJECT_ROOT / "work" / "jshis_official_hazard_response_uncertainty_map.py"

FEATURE_SETS: dict[str, list[str]] = {
    "d1400_avs30": ["log_d1400", "log_avs30"],
    "vs_only": ["log_vs10", "log_vs20", "log_avs30"],
    "depth_only": ["log_d1100", "log_d1400", "log_d1700", "log_d2100", "log_dbase"],
    "location_only": ["lon", "lat"],
    "site_no_location": [
        "log_vs10",
        "log_vs20",
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
    "d1400_avs30": [],
    "vs_only": [],
    "depth_only": [],
    "location_only": [],
    "site_no_location": ["network_label", "installation_situation_id"],
    "full_site_location": ["network_label", "installation_situation_id"],
}

PERMUTATION_GROUPS: dict[str, list[str]] = {
    "shallow_velocity": ["log_vs10", "log_vs20", "log_avs30"],
    "basin_depth": ["log_d1100", "log_d1400", "log_d1700", "log_d2100", "log_dbase"],
    "location": ["lon", "lat"],
    "volcanic_front_distance": ["dist_vf_mf13_nejapan", "dist_vf_mf13_swjapan"],
    "station_metadata": ["elevation", "sensor_depth_glminus", "network_label", "installation_situation_id"],
}

TARGET_SHORT = {
    "PGA RotD50": "PGA",
    "SA(0.3s) RotD50": "SA(0.3 s)",
    "SA(1.0s) RotD50": "SA(1.0 s)",
    "SA(3.0s) RotD50": "SA(3.0 s)",
}


def gradient_boosted(seed: int) -> HistGradientBoostingRegressor:
    return HistGradientBoostingRegressor(
        loss="squared_error",
        learning_rate=0.04,
        max_iter=140,
        max_leaf_nodes=11,
        min_samples_leaf=35,
        l2_regularization=0.1,
        random_state=seed,
    )


def build_custom_pipeline(feature_set: str, seed: int):
    # Reuse the base pipeline builder by temporarily mapping feature-set names.
    # The local preprocessing follows the same imputation/scaling logic.
    from sklearn.compose import ColumnTransformer
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler

    numeric = FEATURE_SETS[feature_set]
    categorical = CATEGORICAL_SETS[feature_set]
    transformers = [
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
    pre = ColumnTransformer(transformers=transformers)
    return Pipeline([("preprocess", pre), ("model", gradient_boosted(seed))])


def split_indices(frame: pd.DataFrame, split: str, seed: int) -> list[tuple[int, np.ndarray, np.ndarray]]:
    if split == "station_random":
        kf = KFold(n_splits=5, shuffle=True, random_state=seed)
        return [(fold, train, test) for fold, (train, test) in enumerate(kf.split(frame))]
    folds = make_spatial_folds(frame, n_splits=5, seed=seed)
    return [(int(fold), np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)) for fold in sorted(np.unique(folds))]


def run_ablation(frame: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    pred_rows: list[pd.DataFrame] = []
    for target in TARGET_ORDER:
        target_frame = frame[frame["target_label"].astype(str).eq(target)].copy().reset_index(drop=True)
        y = target_frame["mean_residual"].to_numpy(dtype=float)
        weights = target_frame["n_records"].to_numpy(dtype=float)
        for split in ["spatial_block"]:
            for fold, train_idx, test_idx in split_indices(target_frame, split, seed):
                y_test = y[test_idx]
                w_test = weights[test_idx]
                baseline_rmse = weighted_rmse(y_test, np.zeros_like(y_test), w_test)
                baseline_mae = weighted_mae(y_test, np.zeros_like(y_test), w_test)
                fold_pred = target_frame.iloc[test_idx][
                    ["target", "target_label", "period", "siteid2", "site_code", "network_label", "n_records", "mean_residual"]
                ].copy()
                fold_pred["split"] = split
                fold_pred["fold"] = fold
                for feature_set in FEATURE_SETS:
                    pipe = build_custom_pipeline(feature_set, seed)
                    pipe.fit(target_frame.iloc[train_idx], y[train_idx], model__sample_weight=weights[train_idx])
                    pred = pipe.predict(target_frame.iloc[test_idx])
                    rmse = weighted_rmse(y_test, pred, w_test)
                    mae = weighted_mae(y_test, pred, w_test)
                    rows.append(
                        {
                            "target_label": target,
                            "split": split,
                            "fold": fold,
                            "feature_set": feature_set,
                            "n_train_sites": len(train_idx),
                            "n_test_sites": len(test_idx),
                            "weighted_mae": mae,
                            "weighted_rmse": rmse,
                            "mae_reduction_pct": 100.0 * (baseline_mae - mae) / baseline_mae,
                            "rmse_reduction_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse,
                            "weighted_r2_vs_zero": 1.0 - (rmse**2 / baseline_rmse**2),
                        }
                    )
                    if feature_set == "full_site_location":
                        fold_pred["full_site_location_prediction"] = pred
                        fold_pred["remaining_residual_after_prediction"] = y_test - pred
                pred_rows.append(fold_pred)
    fold_metrics = pd.DataFrame(rows)
    summary = (
        fold_metrics.groupby(["target_label", "split", "feature_set"], observed=True)
        .agg(
            folds=("fold", "nunique"),
            n_test_sites_mean=("n_test_sites", "mean"),
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
    summary["target_label"] = pd.Categorical(summary["target_label"], TARGET_ORDER, ordered=True)
    summary = summary.sort_values(["split", "target_label", "feature_set"]).reset_index(drop=True)
    predictions = pd.concat(pred_rows, ignore_index=True)
    return summary, fold_metrics, predictions


def run_grouped_permutation_importance(frame: pd.DataFrame, seed: int, repeats: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for target in TARGET_ORDER:
        target_frame = frame[frame["target_label"].astype(str).eq(target)].copy().reset_index(drop=True)
        y = target_frame["mean_residual"].to_numpy(dtype=float)
        weights = target_frame["n_records"].to_numpy(dtype=float)
        for fold, train_idx, test_idx in split_indices(target_frame, "spatial_block", seed):
            pipe = build_custom_pipeline("full_site_location", seed)
            pipe.fit(target_frame.iloc[train_idx], y[train_idx], model__sample_weight=weights[train_idx])
            test_frame = target_frame.iloc[test_idx].copy()
            y_test = y[test_idx]
            w_test = weights[test_idx]
            baseline_pred = pipe.predict(test_frame)
            baseline_rmse = weighted_rmse(y_test, baseline_pred, w_test)
            for group, cols in PERMUTATION_GROUPS.items():
                increases = []
                for _ in range(repeats):
                    permuted = test_frame.copy()
                    shuffled_idx = rng.permutation(permuted.index.to_numpy())
                    for col in cols:
                        permuted[col] = permuted.loc[shuffled_idx, col].to_numpy()
                    pred = pipe.predict(permuted)
                    perm_rmse = weighted_rmse(y_test, pred, w_test)
                    increases.append(100.0 * (perm_rmse - baseline_rmse) / baseline_rmse)
                rows.append(
                    {
                        "target_label": target,
                        "fold": fold,
                        "group": group,
                        "baseline_weighted_rmse": baseline_rmse,
                        "rmse_increase_pct_mean": float(np.mean(increases)),
                        "rmse_increase_pct_std": float(np.std(increases, ddof=1)) if len(increases) > 1 else 0.0,
                    }
                )
    result = pd.DataFrame(rows)
    summary = (
        result.groupby(["target_label", "group"], observed=True)
        .agg(
            folds=("fold", "nunique"),
            rmse_increase_pct_mean=("rmse_increase_pct_mean", "mean"),
            rmse_increase_pct_std=("rmse_increase_pct_mean", "std"),
        )
        .reset_index()
    )
    summary["target_label"] = pd.Categorical(summary["target_label"], TARGET_ORDER, ordered=True)
    return summary.sort_values(["target_label", "group"]).reset_index(drop=True)


def run_network_transfer(frame: pd.DataFrame, seed: int) -> pd.DataFrame:
    rows: list[dict] = []
    for target in TARGET_ORDER:
        target_frame = frame[frame["target_label"].astype(str).eq(target)].copy().reset_index(drop=True)
        for train_network, test_network in [("K-NET", "KiK-net"), ("KiK-net", "K-NET")]:
            train_idx = target_frame.index[target_frame["network_label"].eq(train_network)].to_numpy()
            test_idx = target_frame.index[target_frame["network_label"].eq(test_network)].to_numpy()
            if len(train_idx) == 0 or len(test_idx) == 0:
                continue
            y = target_frame["mean_residual"].to_numpy(dtype=float)
            weights = target_frame["n_records"].to_numpy(dtype=float)
            pipe = build_custom_pipeline("full_site_location", seed)
            pipe.fit(target_frame.iloc[train_idx], y[train_idx], model__sample_weight=weights[train_idx])
            pred = pipe.predict(target_frame.iloc[test_idx])
            y_test = y[test_idx]
            w_test = weights[test_idx]
            baseline_rmse = weighted_rmse(y_test, np.zeros_like(y_test), w_test)
            baseline_mae = weighted_mae(y_test, np.zeros_like(y_test), w_test)
            rmse = weighted_rmse(y_test, pred, w_test)
            mae = weighted_mae(y_test, pred, w_test)
            rows.append(
                {
                    "target_label": target,
                    "train_network": train_network,
                    "test_network": test_network,
                    "n_train_sites": len(train_idx),
                    "n_test_sites": len(test_idx),
                    "weighted_mae": mae,
                    "weighted_rmse": rmse,
                    "mae_reduction_pct": 100.0 * (baseline_mae - mae) / baseline_mae,
                    "rmse_reduction_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse,
                    "weighted_r2_vs_zero": 1.0 - (rmse**2 / baseline_rmse**2),
                }
            )
    result = pd.DataFrame(rows)
    result["target_label"] = pd.Categorical(result["target_label"], TARGET_ORDER, ordered=True)
    return result.sort_values(["target_label", "train_network"]).reset_index(drop=True)


def build_engineering_examples(predictions: pd.DataFrame) -> pd.DataFrame:
    sa3 = predictions[
        predictions["split"].eq("spatial_block") & predictions["target_label"].astype(str).eq("SA(3.0s) RotD50")
    ].copy()
    sa3 = sa3.dropna(subset=["full_site_location_prediction"])
    sa3["predicted_factor"] = np.power(10.0, sa3["full_site_location_prediction"])
    sa3["abs_prediction"] = sa3["full_site_location_prediction"].abs()
    high = sa3.nlargest(6, "full_site_location_prediction")
    low = sa3.nsmallest(6, "full_site_location_prediction")
    out = pd.concat([high, low], ignore_index=True)
    out = out[
        [
            "siteid2",
            "site_code",
            "network_label",
            "n_records",
            "mean_residual",
            "full_site_location_prediction",
            "predicted_factor",
            "remaining_residual_after_prediction",
        ]
    ].copy()
    out = out.rename(
        columns={
            "mean_residual": "heldout_station_residual_log10",
            "full_site_location_prediction": "predicted_station_correction_log10",
            "remaining_residual_after_prediction": "residual_after_correction_log10",
        }
    )
    return out


def write_report(
    ablation_summary: pd.DataFrame,
    importance: pd.DataFrame,
    transfer: pd.DataFrame,
    examples: pd.DataFrame,
    output_path: Path,
) -> None:
    def md(frame: pd.DataFrame, cols: list[str], floatfmt: str = ".2f") -> str:
        sub = frame[cols].copy()
        lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, row in sub.iterrows():
            vals = []
            for col in cols:
                value = row[col]
                if isinstance(value, (float, np.floating)):
                    vals.append(format(float(value), floatfmt))
                else:
                    vals.append(str(value))
            lines.append("| " + " | ".join(vals) + " |")
        return "\n".join(lines)

    spatial = ablation_summary[ablation_summary["split"].eq("spatial_block")]
    key = spatial[
        spatial["target_label"].astype(str).isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])
        & spatial["feature_set"].isin(["d1400_avs30", "vs_only", "depth_only", "location_only", "full_site_location"])
    ].copy()
    transfer_key = transfer[transfer["target_label"].astype(str).isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])]
    lines = [
        "# Station Residual Model Extensions",
        "",
        "## Ablation summary",
        "",
        md(
            key,
            [
                "target_label",
                "feature_set",
                "rmse_reduction_pct_mean",
                "rmse_reduction_pct_std",
                "mae_reduction_pct_mean",
            ],
        ),
        "",
        "## Grouped permutation importance",
        "",
        md(
            importance.sort_values(["target_label", "rmse_increase_pct_mean"], ascending=[True, False]),
            ["target_label", "group", "rmse_increase_pct_mean", "rmse_increase_pct_std"],
        ),
        "",
        "## K-NET/KiK-net transfer",
        "",
        md(
            transfer_key,
            ["target_label", "train_network", "test_network", "n_train_sites", "n_test_sites", "rmse_reduction_pct"],
        ),
        "",
        "## Engineering-scale examples for SA(3.0 s)",
        "",
        md(
            examples,
            [
                "site_code",
                "network_label",
                "n_records",
                "heldout_station_residual_log10",
                "predicted_station_correction_log10",
                "predicted_factor",
                "residual_after_correction_log10",
            ],
            floatfmt=".3f",
        ),
        "",
        "## Interpretation",
        "",
        "The ablation separates basin-depth, shallow-velocity, location, and full site-plus-location information. The grouped permutation analysis evaluates which predictor groups carry the held-out predictive signal. The network-transfer test is stricter than random station splits because training and testing use different Japanese strong-motion networks. The engineering examples translate log10 station corrections into multiplicative spectral-amplitude factors.",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_extension_figure(
    ablation: pd.DataFrame,
    importance: pd.DataFrame,
    transfer: pd.DataFrame,
    examples: pd.DataFrame,
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.8, 8.8))
    grid = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.32)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    # A: spatial-block ablation for long-period targets
    feature_order = ["d1400_avs30", "vs_only", "depth_only", "location_only", "full_site_location"]
    feature_labels = ["D1400+AVS30", "VS only", "Depth only", "Location only", "Full"]
    colors = {"SA(1.0s) RotD50": "#F58518", "SA(3.0s) RotD50": "#E45756"}
    spatial = ablation[
        ablation["split"].eq("spatial_block")
        & ablation["target_label"].astype(str).isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])
    ].copy()
    x = np.arange(len(feature_order))
    width = 0.36
    for offset, target in [(-width / 2, "SA(1.0s) RotD50"), (width / 2, "SA(3.0s) RotD50")]:
        sub = spatial[spatial["target_label"].astype(str).eq(target)].set_index("feature_set").loc[feature_order]
        ax_a.bar(x + offset, sub["rmse_reduction_pct_mean"], width, color=colors[target], label=TARGET_SHORT[target])
    ax_a.axhline(0, color="black", lw=0.8)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(feature_labels, rotation=25, ha="right")
    ax_a.set_ylabel("Spatial-block RMSE reduction (%)")
    ax_a.set_title("A Feature-set ablation")
    ax_a.legend(frameon=False)
    ax_a.grid(axis="y", color="#e6e6e6")

    # B: grouped permutation importance heatmap
    group_order = ["basin_depth", "shallow_velocity", "location", "volcanic_front_distance", "station_metadata"]
    imp = importance.copy()
    pivot = imp.pivot(index="group", columns="target_label", values="rmse_increase_pct_mean").loc[group_order, TARGET_ORDER]
    im = ax_b.imshow(pivot.to_numpy(), cmap="YlOrRd", aspect="auto")
    ax_b.set_xticks(np.arange(len(TARGET_ORDER)))
    ax_b.set_xticklabels([TARGET_SHORT[t] for t in TARGET_ORDER], rotation=25, ha="right")
    ax_b.set_yticks(np.arange(len(group_order)))
    ax_b.set_yticklabels([g.replace("_", " ") for g in group_order])
    ax_b.set_title("B Grouped permutation importance")
    for i, group in enumerate(group_order):
        for j, target in enumerate(TARGET_ORDER):
            ax_b.text(j, i, f"{pivot.loc[group, target]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.04, label="RMSE increase (%)")

    # C: network transfer
    transfer_targets = transfer[transfer["target_label"].astype(str).isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])].copy()
    labels = []
    vals = []
    bar_colors = []
    for _, row in transfer_targets.iterrows():
        labels.append(f"{TARGET_SHORT[str(row['target_label'])]}\n{row['train_network']}→{row['test_network']}")
        vals.append(row["rmse_reduction_pct"])
        bar_colors.append(colors[str(row["target_label"])])
    ax_c.bar(np.arange(len(vals)), vals, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax_c.axhline(0, color="black", lw=0.8)
    ax_c.set_xticks(np.arange(len(vals)))
    ax_c.set_xticklabels(labels, rotation=25, ha="right")
    ax_c.set_ylabel("Network-transfer RMSE reduction (%)")
    ax_c.set_title("C K-NET/KiK-net transfer")
    ax_c.grid(axis="y", color="#e6e6e6")

    # D: engineering examples
    ex = examples.copy()
    ex["label"] = ex["site_code"] + "\n" + ex["network_label"]
    ex = pd.concat([ex.nlargest(4, "predicted_station_correction_log10"), ex.nsmallest(4, "predicted_station_correction_log10")])
    vals = ex["predicted_factor"].to_numpy()
    ax_d.bar(np.arange(len(ex)), vals, color=np.where(vals >= 1.0, "#E45756", "#4C78A8"), edgecolor="black", linewidth=0.5)
    ax_d.axhline(1.0, color="black", lw=0.8)
    ax_d.set_xticks(np.arange(len(ex)))
    ax_d.set_xticklabels(ex["label"], rotation=35, ha="right")
    ax_d.set_ylabel("Predicted SA(3.0 s) factor")
    ax_d.set_title("D Engineering-scale station corrections")
    ax_d.grid(axis="y", color="#e6e6e6")

    fig.suptitle("Ablation, interpretability, transfer, and engineering scale of the station-residual model", fontsize=13)
    fig.savefig(figure_dir / "cee_fig7_station_model_extensions.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / "cee_fig7_station_model_extensions.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_repro_pack(output_dir: Path, figure_dir: Path, pack_dir: Path) -> None:
    if pack_dir.exists():
        shutil.rmtree(pack_dir)
    (pack_dir / "scripts").mkdir(parents=True)
    (pack_dir / "outputs").mkdir()
    (pack_dir / "figures").mkdir()
    shutil.copy2(BASE_SCRIPT, pack_dir / "scripts" / BASE_SCRIPT.name)
    shutil.copy2(SCRIPT_PATH, pack_dir / "scripts" / SCRIPT_PATH.name)
    if HAZARD_SCRIPT.exists():
        shutil.copy2(HAZARD_SCRIPT, pack_dir / "scripts" / HAZARD_SCRIPT.name)
    if OFFICIAL_HAZARD_SCRIPT.exists():
        shutil.copy2(OFFICIAL_HAZARD_SCRIPT, pack_dir / "scripts" / OFFICIAL_HAZARD_SCRIPT.name)
    if OFFICIAL_UNCERTAINTY_SCRIPT.exists():
        shutil.copy2(OFFICIAL_UNCERTAINTY_SCRIPT, pack_dir / "scripts" / OFFICIAL_UNCERTAINTY_SCRIPT.name)
    for name in [
        "jshis_nonergodic_station_model.md",
        "jshis_nonergodic_station_model_summary.csv",
        "jshis_nonergodic_station_model_fold_metrics.csv",
        "jshis_station_model_ablation_summary.csv",
        "jshis_station_model_ablation_fold_metrics.csv",
        "jshis_station_model_grouped_permutation_importance.csv",
        "jshis_station_model_network_transfer.csv",
        "jshis_station_engineering_examples.csv",
        "jshis_station_model_extensions.md",
        "jshis_station_hazard_impact_by_station.csv",
        "jshis_station_hazard_impact_summary.csv",
        "jshis_station_hazard_impact_examples.csv",
        "jshis_station_hazard_curve_points.csv",
        "jshis_station_hazard_impact.md",
        "jshis_official_hazard_response_examples.csv",
        "jshis_official_hazard_curve_api_examples.csv",
        "jshis_official_response_sa3_station_values.csv",
        "jshis_official_response_uhs_examples.csv",
        "jshis_official_hazard_response_summary.csv",
        "jshis_official_hazard_response_check.md",
        "jshis_official_hazard_response_fold_summary.csv",
        "jshis_official_hazard_response_bootstrap_ci.csv",
        "jshis_official_hazard_response_uncertainty_spatial_check.md",
    ]:
        src = output_dir / name
        if src.exists():
            shutil.copy2(src, pack_dir / "outputs" / name)
    for name in [
        "cee_fig6_nonergodic_station_model.pdf",
        "cee_fig6_nonergodic_station_model.png",
        "cee_fig7_station_model_extensions.pdf",
        "cee_fig7_station_model_extensions.png",
        "cee_fig8_hazard_impact.pdf",
        "cee_fig8_hazard_impact.png",
        "cee_fig9_spatial_uncertainty.pdf",
        "cee_fig9_spatial_uncertainty.png",
    ]:
        src = figure_dir / name
        if src.exists():
            shutil.copy2(src, pack_dir / "figures" / name)
    env = PROJECT_ROOT / "outputs" / "cee_environment_snapshot.yml"
    if env.exists():
        shutil.copy2(env, pack_dir / "environment.yml")
    readme = """# Non-ergodic Station Residual Model Reproducibility Pack

This pack contains scripts and aggregate outputs for the station-held-out residual correction experiments.

## Inputs

The scripts expect the following local aggregate inputs in the parent project:

- `outputs/jshis_mf2013_official_station_residual_terms.csv`
- `outputs/jshis_site_schema_v2024_sub1.csv`

Raw J-SHIS/NIED flatfile records are not redistributed here. Follow the J-SHIS/NIED data-use terms and regenerate the aggregate inputs from the official flatfile when preparing a public release.

## Commands

```bash
conda run -n zhy python work/jshis_nonergodic_station_residual_model.py
conda run -n zhy python work/jshis_nonergodic_station_model_extensions.py
conda run -n zhy python work/jshis_station_hazard_impact.py
conda run -n zhy python work/jshis_official_hazard_response_check.py
conda run -n zhy python work/jshis_official_hazard_response_uncertainty_map.py
```

## Included outputs

- held-out station residual model metrics;
- feature-set ablation metrics;
- grouped permutation importance metrics;
- K-NET/KiK-net transfer metrics;
- engineering-scale station correction examples;
- SA(3.0 s) station-multiplier tables for engineering-scale checks;
- official J-SHIS hazard-curve API and response-spectrum map check for SA(3.0 s);
- spatial-block and bootstrap uncertainty summaries for the official response-spectrum check;
- Figure 6, Figure 7, Figure 8 and Figure 9 source graphics.
"""
    (pack_dir / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--station-terms", type=Path, default=DEFAULT_STATION_TERMS)
    parser.add_argument("--site-schema", type=Path, default=DEFAULT_SITE_SCHEMA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--pack-dir", type=Path, default=PROJECT_ROOT / "outputs" / "nonergodic_station_model_repro_pack_v0_1")
    parser.add_argument("--min-records", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260612)
    parser.add_argument("--importance-repeats", type=int, default=3)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    frame = load_station_dataset(args.station_terms, args.site_schema, args.min_records)

    ablation_summary, ablation_folds, ablation_predictions = run_ablation(frame, args.seed)
    importance = run_grouped_permutation_importance(frame, args.seed, args.importance_repeats)
    transfer = run_network_transfer(frame, args.seed)
    examples = build_engineering_examples(ablation_predictions)

    ablation_summary.to_csv(args.output_dir / "jshis_station_model_ablation_summary.csv", index=False)
    ablation_folds.to_csv(args.output_dir / "jshis_station_model_ablation_fold_metrics.csv", index=False)
    ablation_predictions.to_csv(args.output_dir / "jshis_station_model_ablation_predictions.csv", index=False)
    importance.to_csv(args.output_dir / "jshis_station_model_grouped_permutation_importance.csv", index=False)
    transfer.to_csv(args.output_dir / "jshis_station_model_network_transfer.csv", index=False)
    examples.to_csv(args.output_dir / "jshis_station_engineering_examples.csv", index=False)

    write_report(ablation_summary, importance, transfer, examples, args.output_dir / "jshis_station_model_extensions.md")
    build_extension_figure(ablation_summary, importance, transfer, examples, args.figure_dir)
    build_repro_pack(args.output_dir, args.figure_dir, args.pack_dir)

    print(f"Wrote {args.output_dir / 'jshis_station_model_ablation_summary.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_model_grouped_permutation_importance.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_model_network_transfer.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_engineering_examples.csv'}")
    print(f"Wrote {args.figure_dir / 'cee_fig7_station_model_extensions.pdf'}")
    print(f"Wrote {args.pack_dir}")


if __name__ == "__main__":
    main()
