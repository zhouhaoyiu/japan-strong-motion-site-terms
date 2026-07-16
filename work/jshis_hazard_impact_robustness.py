#!/usr/bin/env python3
"""Quantify uncertainty and independent-model stability of spectrum corrections."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core
from jshis_station_uncertainty_propagation import weighted_quantile


PRIMARY_PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv"
CROSS_NETWORK_PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_cross_network_transfer_predictions.csv"

OUT_INTERVAL_CURVE = core.SUPPLEMENT_DIR / "jshis_hazard_impact_interval_robustness.csv"
OUT_MODEL_STATIONS = core.SUPPLEMENT_DIR / "jshis_hazard_impact_independent_model_stations.csv"
OUT_MODEL_SUMMARY = core.SUPPLEMENT_DIR / "jshis_hazard_impact_independent_model_summary.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_hazard_impact_robustness.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "supplementary_figure_hazard_impact_robustness.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "supplementary_figure_hazard_impact_robustness.png"


def primary_frame(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame = frame[frame["model"].eq("physical_spatial_hgb")].copy()
    if "residual_model" in frame.columns:
        frame = frame[frame["residual_model"].eq("mf2013_site")].copy()
    if frame.duplicated(["siteid2", "period_s"]).any():
        raise ValueError("Primary prediction table contains duplicate station-period rows")
    return frame


def interval_robustness(
    predictions: pd.DataFrame,
    confidence_levels: list[float],
    material_factor: float,
) -> pd.DataFrame:
    rows = []
    material_log10 = float(np.log10(material_factor))
    for confidence in confidence_levels:
        lower_quantile = (1.0 - confidence) / 2.0
        upper_quantile = 1.0 - lower_quantile
        for period_s, period_frame in predictions.groupby("period_s", sort=True):
            period_frame = period_frame.copy().reset_index(drop=True)
            observed = period_frame["station_effect_log10"].to_numpy(float)
            point = period_frame["centered_oof_prediction_log10"].to_numpy(float)
            weights = period_frame["n_records"].to_numpy(float)
            folds = period_frame["fold"].to_numpy(int)
            lower = np.full(len(period_frame), np.nan)
            upper = np.full(len(period_frame), np.nan)
            for fold in sorted(np.unique(folds)):
                calibration = folds != fold
                held = folds == fold
                errors = observed[calibration] - point[calibration]
                q_low, q_high = weighted_quantile(
                    errors,
                    [lower_quantile, upper_quantile],
                    weights[calibration],
                )
                lower[held] = point[held] + q_low
                upper[held] = point[held] + q_high
            covered = (observed >= lower) & (observed <= upper)
            robust_increase = lower > 0
            robust_decrease = upper < 0
            material_increase = lower > material_log10
            material_decrease = upper < -material_log10
            robust = robust_increase | robust_decrease
            material = material_increase | material_decrease
            rows.append(
                {
                    "period_s": period_s,
                    "confidence_level_pct": 100.0 * confidence,
                    "n_stations": len(period_frame),
                    "empirical_coverage_pct": 100.0 * covered.mean(),
                    "record_weighted_coverage_pct": 100.0
                    * np.average(covered, weights=weights),
                    "robust_direction_fraction_pct": 100.0 * robust.mean(),
                    "record_weighted_robust_direction_fraction_pct": 100.0
                    * np.average(robust, weights=weights),
                    "robust_increase_fraction_pct": 100.0 * robust_increase.mean(),
                    "robust_decrease_fraction_pct": 100.0 * robust_decrease.mean(),
                    "robust_material_fraction_pct": 100.0 * material.mean(),
                    "record_weighted_robust_material_fraction_pct": 100.0
                    * np.average(material, weights=weights),
                    "material_factor": material_factor,
                }
            )
    return pd.DataFrame(rows)


def aggregate_cross_network(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame = frame[
        frame["source_network"].eq("K-NET")
        & frame["target_network"].eq("KiK-net")
        & frame["model"].eq("physical_spatial_hgb")
    ].copy()
    rows = []
    identity = ["period_s", "period_code", "siteid2", "site_code", "lon", "lat", "avs30", "d1400"]
    for keys, block in frame.groupby(identity, dropna=False, sort=False):
        weights = block["n_test_records"].to_numpy(float)
        row = dict(zip(identity, keys, strict=True))
        row["cross_network_prediction_log10"] = float(
            np.average(block["frozen_prediction_log10"], weights=weights)
        )
        row["cross_network_test_records"] = int(weights.sum())
        rows.append(row)
    return pd.DataFrame(rows)


def independent_model_agreement(
    primary: pd.DataFrame,
    cross_network: pd.DataFrame,
    material_factor: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    keep = [
        "period_s",
        "siteid2",
        "n_records",
        "centered_oof_prediction_log10",
        "station_effect_log10",
    ]
    target = primary[primary["network_label"].eq("KiK-net")][keep].copy()
    target = target.rename(
        columns={
            "n_records": "primary_record_weight",
            "centered_oof_prediction_log10": "primary_prediction_log10",
            "station_effect_log10": "observed_station_term_log10",
        }
    )
    stations = cross_network.merge(
        target,
        on=["period_s", "siteid2"],
        how="inner",
        validate="one_to_one",
    )
    stations["primary_multiplier"] = 10.0 ** stations["primary_prediction_log10"]
    stations["cross_network_multiplier"] = 10.0 ** stations["cross_network_prediction_log10"]
    stations["absolute_model_difference_log10"] = (
        stations["primary_prediction_log10"] - stations["cross_network_prediction_log10"]
    ).abs()
    stations["multiplier_ratio_cross_over_primary"] = (
        stations["cross_network_multiplier"] / stations["primary_multiplier"]
    )
    stations["direction_agreement"] = np.sign(stations["primary_prediction_log10"]).eq(
        np.sign(stations["cross_network_prediction_log10"])
    )
    threshold = float(np.log10(material_factor))
    stations["material_increase_agreement"] = stations["primary_prediction_log10"].gt(
        threshold
    ) & stations["cross_network_prediction_log10"].gt(threshold)
    stations["material_decrease_agreement"] = stations["primary_prediction_log10"].lt(
        -threshold
    ) & stations["cross_network_prediction_log10"].lt(-threshold)
    stations["material_direction_agreement"] = (
        stations["material_increase_agreement"]
        | stations["material_decrease_agreement"]
    )

    rows = []
    for period_s, block in stations.groupby("period_s", sort=True):
        weights = block["primary_record_weight"].to_numpy(float)
        primary_values = block["primary_prediction_log10"].to_numpy(float)
        cross_values = block["cross_network_prediction_log10"].to_numpy(float)
        rows.append(
            {
                "period_s": period_s,
                "n_stations": len(block),
                "independent_model_pearson": block["primary_prediction_log10"].corr(
                    block["cross_network_prediction_log10"]
                ),
                "record_weighted_model_pearson": float(
                    np.cov(primary_values, cross_values, aweights=weights)[0, 1]
                    / np.sqrt(
                        np.cov(primary_values, cross_values, aweights=weights)[0, 0]
                        * np.cov(primary_values, cross_values, aweights=weights)[1, 1]
                    )
                ),
                "direction_agreement_pct": 100.0 * block["direction_agreement"].mean(),
                "record_weighted_direction_agreement_pct": 100.0
                * np.average(block["direction_agreement"], weights=weights),
                "material_direction_agreement_pct": 100.0
                * block["material_direction_agreement"].mean(),
                "record_weighted_material_direction_agreement_pct": 100.0
                * np.average(block["material_direction_agreement"], weights=weights),
                "median_absolute_model_difference_log10": block[
                    "absolute_model_difference_log10"
                ].median(),
                "q95_absolute_model_difference_log10": block[
                    "absolute_model_difference_log10"
                ].quantile(0.95),
                "multiplier_ratio_q05": block["multiplier_ratio_cross_over_primary"].quantile(
                    0.05
                ),
                "multiplier_ratio_q50": block["multiplier_ratio_cross_over_primary"].median(),
                "multiplier_ratio_q95": block["multiplier_ratio_cross_over_primary"].quantile(
                    0.95
                ),
                "material_factor": material_factor,
            }
        )
    return stations, pd.DataFrame(rows)


def plot_results(
    interval_curve: pd.DataFrame,
    stations: pd.DataFrame,
    model_summary: pd.DataFrame,
    material_factor: float,
) -> None:
    selected_periods = [0.5, 1.0, 2.0, 3.0, 5.0]
    colors = {
        0.5: "#CC79A7",
        1.0: "#E69F00",
        2.0: "#009E73",
        3.0: "#0072B2",
        5.0: "#D55E00",
    }
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot([45, 100], [45, 100], color="#555555", ls="--", lw=1.0, label="Nominal")
    for period_s in selected_periods:
        block = interval_curve[interval_curve["period_s"].eq(period_s)].sort_values(
            "confidence_level_pct"
        )
        ax.plot(
            block["confidence_level_pct"],
            block["empirical_coverage_pct"],
            "o-",
            color=colors[period_s],
            label=f"{period_s:g} s",
        )
    ax.set(
        xlabel="Nominal interval coverage (%)",
        ylabel="Empirical held-block coverage (%)",
        title="a  Spatial-block interval calibration",
    )
    ax.legend(frameon=False, ncol=2, fontsize=8)

    ax = axes[0, 1]
    for period_s in selected_periods:
        block = interval_curve[interval_curve["period_s"].eq(period_s)].sort_values(
            "confidence_level_pct"
        )
        ax.plot(
            block["confidence_level_pct"],
            block["robust_direction_fraction_pct"],
            "o-",
            color=colors[period_s],
            label=f"{period_s:g} s",
        )
    ax.set(
        xlabel="Prediction interval coverage (%)",
        ylabel="Stations with interval excluding unity (%)",
        title="b  Directional robustness under prediction error",
    )

    ax = axes[1, 0]
    block3 = stations[stations["period_s"].eq(3.0)].copy()
    categories = [
        (block3["material_decrease_agreement"], "Both decrease", "#2166AC"),
        (block3["material_increase_agreement"], "Both increase", "#B2182B"),
        (~block3["material_direction_agreement"] & block3["direction_agreement"], "Same sign, <10% in one model", "#BDBDBD"),
        (~block3["direction_agreement"], "Direction differs", "#FDD49E"),
    ]
    for mask, label, color in categories:
        ax.scatter(
            block3.loc[mask, "lon"],
            block3.loc[mask, "lat"],
            s=13,
            color=color,
            linewidths=0,
            label=label,
        )
    ax.set(
        xlabel="Longitude (°E)",
        ylabel="Latitude (°N)",
        title=f"c  Independent-model agreement at 3.0 s (factor {material_factor:.1f})",
    )
    ax.legend(frameon=False, fontsize=7, loc="lower right")

    ax = axes[1, 1]
    ax.scatter(
        block3["primary_prediction_log10"],
        block3["cross_network_prediction_log10"],
        s=np.clip(block3["primary_record_weight"], 8, 70),
        alpha=0.5,
        color="#2878B5",
        linewidths=0,
    )
    limits = np.quantile(
        np.r_[block3["primary_prediction_log10"], block3["cross_network_prediction_log10"]],
        [0.005, 0.995],
    )
    ax.plot(limits, limits, color="#555555", ls="--", lw=1.0)
    ax.text(0.96, 0.04, "Dashed: equality", transform=ax.transAxes, ha="right", color="#555555")
    row3 = model_summary[model_summary["period_s"].eq(3.0)].iloc[0]
    ax.text(
        0.04,
        0.96,
        (
            f"Pearson $r$ = {row3['independent_model_pearson']:.3f}\n"
            f"Direction agreement = {row3['direction_agreement_pct']:.1f}%\n"
            f"Both exceed 10% = {row3['material_direction_agreement_pct']:.1f}%"
        ),
        transform=ax.transAxes,
        va="top",
    )
    ax.set(
        xlim=limits,
        ylim=limits,
        xlabel="Primary spatial-block prediction (log$_{10}$)",
        ylabel="Frozen K-NET prediction (log$_{10}$)",
        title="d  Independent correction fields at KiK-net sites",
    )

    for axis in axes.flat:
        axis.grid(alpha=0.18, linewidth=0.6)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    primary = primary_frame(args.primary_predictions)
    confidence_levels = [float(value) for value in args.confidence_levels.split(",")]
    interval_curve = interval_robustness(
        primary,
        confidence_levels,
        args.material_factor,
    )
    cross_network = aggregate_cross_network(args.cross_network_predictions)
    stations, model_summary = independent_model_agreement(
        primary,
        cross_network,
        args.material_factor,
    )
    interval_curve.to_csv(OUT_INTERVAL_CURVE, index=False)
    stations.to_csv(OUT_MODEL_STATIONS, index=False)
    model_summary.to_csv(OUT_MODEL_SUMMARY, index=False)
    plot_results(interval_curve, stations, model_summary, args.material_factor)

    interval3 = interval_curve[interval_curve["period_s"].eq(3.0)].set_index(
        "confidence_level_pct"
    )
    model3 = model_summary[model_summary["period_s"].eq(3.0)].iloc[0]
    lines = [
        "# Hazard-impact uncertainty and independent-model stability",
        "",
        "Spatial-block error intervals and independent network-transfer predictions address different uncertainty questions and are reported separately.",
        "",
        f"- SA(3.0 s) empirical coverage at nominal 80% and 90%: {interval3.loc[80.0, 'empirical_coverage_pct']:.1f}% and {interval3.loc[90.0, 'empirical_coverage_pct']:.1f}%.",
        f"- SA(3.0 s) intervals excluding unity at nominal 80% and 90%: {interval3.loc[80.0, 'robust_direction_fraction_pct']:.1f}% and {interval3.loc[90.0, 'robust_direction_fraction_pct']:.1f}% of stations.",
        f"- Independent K-NET-trained and primary correction-field correlation at 3.0 s: {model3['independent_model_pearson']:.3f} across {int(model3['n_stations']):,} KiK-net sites.",
        f"- Independent-model direction agreement: {model3['direction_agreement_pct']:.1f}%.",
        f"- Independent-model agreement on changes exceeding a factor of {args.material_factor:.1f}: {model3['material_direction_agreement_pct']:.1f}%.",
        f"- Median and 95th-percentile absolute difference between the two correction fields: {model3['median_absolute_model_difference_log10']:.3f} and {model3['q95_absolute_model_difference_log10']:.3f} log10 units.",
        "",
        "Prediction intervals include station-model error under spatial blocking. Independent-model agreement measures reproducibility under a disjoint network and earthquake set. Neither quantity includes source-occurrence, ground-motion-backbone, or J-SHIS logic-tree uncertainty.",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT_INTERVAL_CURVE}", flush=True)
    print(f"wrote {OUT_MODEL_STATIONS}", flush=True)
    print(f"wrote {OUT_MODEL_SUMMARY}", flush=True)
    print(f"wrote {OUT_AUDIT}", flush=True)
    print(f"wrote {OUT_FIGURE_PDF}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-predictions", type=Path, default=PRIMARY_PREDICTIONS)
    parser.add_argument(
        "--cross-network-predictions",
        type=Path,
        default=CROSS_NETWORK_PREDICTIONS,
    )
    parser.add_argument("--confidence-levels", default="0.50,0.60,0.70,0.80,0.90,0.95")
    parser.add_argument("--material-factor", type=float, default=1.10)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
