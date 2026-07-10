#!/usr/bin/env python3
"""Predict equal-stratum station fields with the primary spatial validation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core


PATH_TERMS = core.SUPPLEMENT_DIR / "jshis_path_stratification_station_terms.csv"
PRIMARY_TERMS = core.OUT_STATION_TERMS
PRIMARY_METRICS = core.OUT_MODEL_METRICS
OUT_TERMS = core.SUPPLEMENT_DIR / "jshis_balanced_station_terms.csv"
OUT_PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_balanced_station_model_predictions.csv"
OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_balanced_station_model_metrics.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_balanced_station_model_audit.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "supplementary_figure_balanced_station_model.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "supplementary_figure_balanced_station_model.png"

TARGETS = {
    "hypocentral_azimuth": "azimuth_balanced",
    "fault_distance": "distance_balanced",
    "source_type": "source_balanced",
}


def build_balanced_terms(
    path_terms: pd.DataFrame,
    primary_terms: pd.DataFrame,
    min_strata: int,
) -> pd.DataFrame:
    metadata = primary_terms[primary_terms["model"].eq("mf2013_site")].drop_duplicates(
        ["siteid2", "period_s"]
    )
    metadata = metadata.drop(
        columns=["n_records", "station_effect_log10", "model", "response_component"],
    )
    blocks = []
    for dimension, target in TARGETS.items():
        block = path_terms[path_terms["stratification"].eq(dimension)]
        balanced = (
            block.groupby(["siteid2", "period_s"], as_index=False)
            .agg(
                station_effect_log10=("aligned_stratum_station_effect_log10", "mean"),
                n_records=("stratum_n_records", "sum"),
                n_strata=("stratum", "nunique"),
            )
        )
        balanced = balanced[balanced["n_strata"].ge(min_strata)].merge(
            metadata,
            on=["siteid2", "period_s"],
            how="left",
            validate="one_to_one",
        )
        balanced["model"] = target
        balanced["response_component"] = "RotD100"
        balanced["stratification"] = dimension
        blocks.append(balanced)
    return pd.concat(blocks, ignore_index=True)


def save_figure(metrics: pd.DataFrame, primary_metrics: pd.DataFrame) -> None:
    selected_overall = metrics[
        metrics["scope"].eq("overall") & metrics["model"].eq("physical_spatial_hgb")
    ].copy()
    selected_folds = metrics[
        metrics["scope"].eq("fold") & metrics["model"].eq("physical_spatial_hgb")
    ].copy()
    primary_overall = primary_metrics[
        primary_metrics["scope"].eq("overall")
        & primary_metrics["model"].eq("physical_spatial_hgb")
    ].copy()
    primary_folds = primary_metrics[
        primary_metrics["scope"].eq("fold")
        & primary_metrics["model"].eq("physical_spatial_hgb")
    ].copy()
    primary_overall["residual_model"] = "primary_full_sample"
    primary_folds["residual_model"] = "primary_full_sample"
    selected_overall = pd.concat([primary_overall, selected_overall], ignore_index=True)
    selected_folds = pd.concat([primary_folds, selected_folds], ignore_index=True)
    labels = {
        "primary_full_sample": "Full sample",
        "azimuth_balanced": "Azimuth balanced",
        "distance_balanced": "Distance balanced",
        "source_balanced": "Source-class balanced",
    }
    colors = {
        "primary_full_sample": "#111111",
        "azimuth_balanced": "#126782",
        "distance_balanced": "#E76F51",
        "source_balanced": "#2A9D8F",
    }
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.7,
            "savefig.dpi": 300,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), constrained_layout=True)
    for target, block in selected_overall.groupby("residual_model"):
        block = block.sort_values("period_s")
        fold_range = (
            selected_folds[selected_folds["residual_model"].eq(target)]
            .groupby("period_s", as_index=False)
            .agg(
                gain_min=("rmse_reduction_vs_zero_pct", "min"),
                gain_max=("rmse_reduction_vs_zero_pct", "max"),
                correlation_min=("observed_predicted_correlation", "min"),
                correlation_max=("observed_predicted_correlation", "max"),
            )
            .sort_values("period_s")
        )
        axes[0].errorbar(
            fold_range["period_s"],
            block["rmse_reduction_vs_zero_pct"],
            yerr=[
                np.maximum(
                    0.0,
                    block["rmse_reduction_vs_zero_pct"].to_numpy()
                    - fold_range["gain_min"].to_numpy(),
                ),
                np.maximum(
                    0.0,
                    fold_range["gain_max"].to_numpy()
                    - block["rmse_reduction_vs_zero_pct"].to_numpy(),
                ),
            ],
            fmt="none",
            color=colors[target],
            alpha=0.35,
            linewidth=0.7,
            capsize=1.5,
        )
        axes[1].errorbar(
            fold_range["period_s"],
            block["observed_predicted_correlation"],
            yerr=[
                np.maximum(
                    0.0,
                    block["observed_predicted_correlation"].to_numpy()
                    - fold_range["correlation_min"].to_numpy(),
                ),
                np.maximum(
                    0.0,
                    fold_range["correlation_max"].to_numpy()
                    - block["observed_predicted_correlation"].to_numpy(),
                ),
            ],
            fmt="none",
            color=colors[target],
            alpha=0.35,
            linewidth=0.7,
            capsize=1.5,
        )
        axes[0].semilogx(
            block["period_s"],
            block["rmse_reduction_vs_zero_pct"],
            marker="o",
            linewidth=1.1,
            color=colors[target],
            label=labels[target],
        )
        axes[1].semilogx(
            block["period_s"],
            block["observed_predicted_correlation"],
            marker="o",
            linewidth=1.1,
            color=colors[target],
            label=labels[target],
        )
    axes[0].axhline(0.0, color="#777777", linewidth=0.7)
    axes[0].set(title="a  Spatial-block RMSE gain", ylabel="Gain versus zero (%)")
    axes[1].set(title="b  Observed-predicted agreement", ylabel="Pearson correlation")
    for ax in axes:
        ax.set_xlabel("Oscillator period (s)")
        ax.set_xticks([0.1, 0.2, 0.5, 1.0, 2.0, 5.0], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
        ax.grid(alpha=0.2, linewidth=0.5)
    axes[1].legend(frameon=False, loc="best")
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, bbox_inches="tight")
    plt.close(fig)


def write_audit(terms: pd.DataFrame, metrics: pd.DataFrame, min_strata: int) -> None:
    overall = metrics[
        metrics["scope"].eq("overall")
        & metrics["model"].eq("physical_spatial_hgb")
        & metrics["period_s"].isin([1.0, 2.0, 3.0])
    ]
    folds = metrics[
        metrics["scope"].eq("fold")
        & metrics["model"].eq("physical_spatial_hgb")
        & metrics["period_s"].isin([1.0, 2.0, 3.0])
    ]
    lines = [
        "# Equal-stratum station-field sensitivity",
        "",
        "## Design",
        "",
        "- Component-aligned station terms are averaged with equal weight across available hypocentral-bearing sectors, distance bins or source classes.",
        f"- A station must contribute to at least {min_strata} strata. The public-variable model uses the same five-block clustering procedure and hyperparameters as the primary analysis. Target strata receive equal weight; validation errors retain station record-count weights.",
        "",
        "## Results at 1-3 s",
        "",
    ]
    for target, label in [
        ("azimuth_balanced", "Azimuth-balanced field"),
        ("distance_balanced", "Distance-balanced field"),
        ("source_balanced", "Source-class-balanced field"),
    ]:
        block = overall[overall["residual_model"].eq(target)]
        fold_min = folds[folds["residual_model"].eq(target)]["rmse_reduction_vs_zero_pct"].min()
        count = terms[terms["model"].eq(target) & terms["period_s"].eq(3.0)]["siteid2"].nunique()
        lines.append(
            f"- {label}: {count:,} stations at 3 s; RMSE gains span "
            f"{block['rmse_reduction_vs_zero_pct'].min():.1f}%--{block['rmse_reduction_vs_zero_pct'].max():.1f}% and correlations span "
            f"{block['observed_predicted_correlation'].min():.3f}--{block['observed_predicted_correlation'].max():.3f}; "
            f"the minimum held-block RMSE gain is {fold_min:.1f}%."
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "Equal-stratum averaging reduces dominance by the observed path distribution. It does not identify a pure site term, because station-specific path interactions can remain within each stratum and coverage is incomplete. These results are a sensitivity analysis of the persistent station-associated component.",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    core.SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    core.FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    path_terms = pd.read_csv(args.path_terms)
    primary_terms = pd.read_csv(args.primary_terms)
    balanced_terms = build_balanced_terms(path_terms, primary_terms, args.min_strata)
    prediction_blocks = []
    metric_blocks = []
    for target in TARGETS.values():
        predictions, metrics = core.cross_validate_station_models(
            balanced_terms,
            min_records=args.min_records,
            n_splits=args.spatial_folds,
            seed=args.seed,
            target_model=target,
        )
        prediction_blocks.append(predictions)
        metric_blocks.append(metrics)
    predictions = pd.concat(prediction_blocks, ignore_index=True)
    metrics = pd.concat(metric_blocks, ignore_index=True)
    balanced_terms.to_csv(OUT_TERMS, index=False)
    predictions.to_csv(OUT_PREDICTIONS, index=False)
    metrics.to_csv(OUT_METRICS, index=False)
    primary_metrics = pd.read_csv(args.primary_metrics)
    save_figure(metrics, primary_metrics)
    write_audit(balanced_terms, metrics, args.min_strata)
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-terms", type=Path, default=PATH_TERMS)
    parser.add_argument("--primary-terms", type=Path, default=PRIMARY_TERMS)
    parser.add_argument("--primary-metrics", type=Path, default=PRIMARY_METRICS)
    parser.add_argument("--min-strata", type=int, default=2)
    parser.add_argument("--min-records", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
