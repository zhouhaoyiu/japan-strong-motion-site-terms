#!/usr/bin/env python3
"""Calibrate spatial-block prediction intervals and propagate them to spectra.

The intervals quantify uncertainty in the out-of-fold MF2013 station-term
prediction. For every held spatial block, signed prediction errors from the
other four blocks define a weighted empirical 90% interval. The resulting
bounds are then applied multiplicatively to the matched J-SHIS response-spectrum
ordinates. They are not a complete uncertainty model for probabilistic seismic
hazard analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core


SUPPLEMENT_DIR = core.SUPPLEMENT_DIR
FIGURE_DIR = core.FIGURE_DIR

PREDICTIONS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv"
HAZARD_VALUES = SUPPLEMENT_DIR / "jshis_event_adjusted_surface_spectrum_values.csv"
CITY_CASES = SUPPLEMENT_DIR / "jshis_event_adjusted_city_nearest_cases.csv"

OUT_INTERVALS = SUPPLEMENT_DIR / "jshis_station_term_prediction_intervals.csv"
OUT_COVERAGE = SUPPLEMENT_DIR / "jshis_station_term_interval_coverage.csv"
OUT_HAZARD_INTERVALS = SUPPLEMENT_DIR / "jshis_surface_spectrum_prediction_intervals.csv"
OUT_HAZARD_SUMMARY = SUPPLEMENT_DIR / "jshis_surface_spectrum_interval_summary.csv"
OUT_AUDIT = SUPPLEMENT_DIR / "jshis_station_uncertainty_propagation.md"
OUT_FIGURE_PDF = FIGURE_DIR / "supplementary_figure_station_uncertainty.pdf"
OUT_FIGURE_PNG = FIGURE_DIR / "supplementary_figure_station_uncertainty.png"


def weighted_quantile(values: np.ndarray, quantiles: list[float], weights: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        raise ValueError("weighted_quantile received no valid values")
    values = values[valid]
    weights = weights[valid]
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    positions = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
    return np.interp(np.asarray(quantiles), positions, values, left=values[0], right=values[-1])


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(np.asarray(values, dtype=float), weights=np.asarray(weights, dtype=float)))


def calibrate_intervals(
    predictions: pd.DataFrame,
    lower_quantile: float,
    upper_quantile: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = predictions[predictions["model"].eq("physical_spatial_hgb")].copy()
    if "residual_model" in frame.columns:
        frame = frame[frame["residual_model"].eq("mf2013_site")].copy()
    required = [
        "siteid2",
        "site_code",
        "network_label",
        "meshcode250_str",
        "lon",
        "lat",
        "avs30",
        "d1400",
        "n_records",
        "station_effect_log10",
        "period_s",
        "period_code",
        "fold",
        "centered_oof_prediction_log10",
    ]
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"prediction table is missing columns: {missing}")
    if frame.duplicated(["siteid2", "period_s"]).any():
        raise ValueError("prediction table has duplicate station-period rows")

    interval_blocks: list[pd.DataFrame] = []
    metric_rows: list[dict[str, float | int | str]] = []
    for period_s, period_frame in frame.groupby("period_s", sort=True):
        period_frame = period_frame.copy().reset_index(drop=True)
        y = period_frame["station_effect_log10"].to_numpy(float)
        point = period_frame["centered_oof_prediction_log10"].to_numpy(float)
        weights = period_frame["n_records"].to_numpy(float)
        folds = period_frame["fold"].to_numpy(int)
        period_rows: list[pd.DataFrame] = []
        for fold in sorted(np.unique(folds)):
            calibration = folds != fold
            held = folds == fold
            errors = y[calibration] - point[calibration]
            q_low, q_high = weighted_quantile(
                errors,
                [lower_quantile, upper_quantile],
                weights[calibration],
            )
            block = period_frame.loc[held, required].copy()
            block["point_prediction_log10"] = point[held]
            block["calibration_error_q_low_log10"] = q_low
            block["calibration_error_q_high_log10"] = q_high
            block["lower_prediction_log10"] = point[held] + q_low
            block["upper_prediction_log10"] = point[held] + q_high
            block["lower_multiplier"] = 10.0 ** block["lower_prediction_log10"]
            block["point_multiplier"] = 10.0 ** block["point_prediction_log10"]
            block["upper_multiplier"] = 10.0 ** block["upper_prediction_log10"]
            block["interval_factor_span"] = block["upper_multiplier"] / block["lower_multiplier"]
            block["covered"] = (
                block["station_effect_log10"].ge(block["lower_prediction_log10"])
                & block["station_effect_log10"].le(block["upper_prediction_log10"])
            )
            interval_blocks.append(block)
            period_rows.append(block)

            covered = block["covered"].to_numpy(float)
            held_weights = block["n_records"].to_numpy(float)
            metric_rows.append(
                {
                    "period_s": period_s,
                    "period_code": str(block["period_code"].iloc[0]),
                    "scope": "fold",
                    "fold": int(fold),
                    "nominal_coverage_pct": 100.0 * (upper_quantile - lower_quantile),
                    "n_calibration_stations": int(calibration.sum()),
                    "n_test_stations": int(held.sum()),
                    "test_record_weight": int(held_weights.sum()),
                    "station_coverage_pct": 100.0 * float(covered.mean()),
                    "record_weighted_coverage_pct": 100.0 * weighted_mean(covered, held_weights),
                    "median_interval_width_log10": float(
                        np.median(block["upper_prediction_log10"] - block["lower_prediction_log10"])
                    ),
                    "median_interval_factor_span": float(np.median(block["interval_factor_span"])),
                    "calibration_error_q_low_log10": q_low,
                    "calibration_error_q_high_log10": q_high,
                }
            )

        combined = pd.concat(period_rows, ignore_index=True)
        combined_covered = combined["covered"].to_numpy(float)
        combined_weights = combined["n_records"].to_numpy(float)
        metric_rows.append(
            {
                "period_s": period_s,
                "period_code": str(combined["period_code"].iloc[0]),
                "scope": "overall",
                "fold": -1,
                "nominal_coverage_pct": 100.0 * (upper_quantile - lower_quantile),
                "n_calibration_stations": int(len(combined) - combined.groupby("fold").size().mean()),
                "n_test_stations": len(combined),
                "test_record_weight": int(combined_weights.sum()),
                "station_coverage_pct": 100.0 * float(combined_covered.mean()),
                "record_weighted_coverage_pct": 100.0 * weighted_mean(combined_covered, combined_weights),
                "median_interval_width_log10": float(
                    np.median(combined["upper_prediction_log10"] - combined["lower_prediction_log10"])
                ),
                "median_interval_factor_span": float(np.median(combined["interval_factor_span"])),
                "calibration_error_q_low_log10": np.nan,
                "calibration_error_q_high_log10": np.nan,
            }
        )
        print(f"interval calibration: period={period_s:g}s stations={len(combined):,}", flush=True)

    intervals = pd.concat(interval_blocks, ignore_index=True)
    intervals = intervals.sort_values(["period_s", "siteid2"]).reset_index(drop=True)
    coverage = pd.DataFrame(metric_rows).sort_values(["period_s", "scope", "fold"]).reset_index(drop=True)
    return intervals, coverage


def propagate_to_spectra(
    hazard_values: pd.DataFrame,
    intervals: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    hazard = hazard_values[hazard_values["model"].eq("physical_spatial_hgb")].copy()
    if "residual_model" in hazard.columns:
        hazard = hazard[hazard["residual_model"].eq("mf2013_site")].copy()
    keep = [
        "siteid2",
        "period_s",
        "point_prediction_log10",
        "lower_prediction_log10",
        "upper_prediction_log10",
        "lower_multiplier",
        "point_multiplier",
        "upper_multiplier",
        "interval_factor_span",
    ]
    result = hazard.merge(intervals[keep], on=["siteid2", "period_s"], how="inner", validate="many_to_one")
    result["adjusted_surface_sa_lower_g"] = result["ergodic_surface_sa_g"] * result["lower_multiplier"]
    result["adjusted_surface_sa_point_g"] = result["ergodic_surface_sa_g"] * result["point_multiplier"]
    result["adjusted_surface_sa_upper_g"] = result["ergodic_surface_sa_g"] * result["upper_multiplier"]
    if not np.allclose(
        result["adjusted_surface_sa_point_g"],
        result["adjusted_surface_sa_g"],
        rtol=1e-9,
        atol=1e-12,
    ):
        raise ValueError("propagated point estimate does not match the primary hazard table")
    if not (
        result["adjusted_surface_sa_lower_g"].le(result["adjusted_surface_sa_point_g"])
        & result["adjusted_surface_sa_point_g"].le(result["adjusted_surface_sa_upper_g"])
    ).all():
        raise ValueError("one or more propagated intervals do not contain the point estimate")

    summary_rows: list[dict[str, float | int | str]] = []
    for (period_s, probability_level), block in result.groupby(["period_s", "probability_level"], sort=True):
        summary_rows.append(
            {
                "period_s": period_s,
                "probability_level": probability_level,
                "n_stations": block["siteid2"].nunique(),
                "median_ergodic_surface_sa_g": float(block["ergodic_surface_sa_g"].median()),
                "median_adjusted_surface_sa_point_g": float(block["adjusted_surface_sa_point_g"].median()),
                "median_adjusted_surface_sa_lower_g": float(block["adjusted_surface_sa_lower_g"].median()),
                "median_adjusted_surface_sa_upper_g": float(block["adjusted_surface_sa_upper_g"].median()),
                "median_point_multiplier": float(block["point_multiplier"].median()),
                "point_multiplier_q05": float(block["point_multiplier"].quantile(0.05)),
                "point_multiplier_q95": float(block["point_multiplier"].quantile(0.95)),
                "median_interval_factor_span": float(block["interval_factor_span"].median()),
            }
        )
    return result, pd.DataFrame(summary_rows)


def plot_figure(
    intervals: pd.DataFrame,
    coverage: pd.DataFrame,
    hazard_intervals: pd.DataFrame,
    city_cases: pd.DataFrame,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.7,
            "savefig.dpi": 300,
        }
    )
    colors = {"point": "#126782", "interval": "#F28E2B", "neutral": "#4D4D4D"}
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2), constrained_layout=True)

    overall = coverage[coverage["scope"].eq("overall")].sort_values("period_s")
    ax = axes[0, 0]
    ax.semilogx(
        overall["period_s"],
        overall["station_coverage_pct"],
        marker="o",
        color=colors["point"],
        label="Stations",
    )
    ax.semilogx(
        overall["period_s"],
        overall["record_weighted_coverage_pct"],
        marker="s",
        color=colors["interval"],
        label="Record weighted",
    )
    ax.axhline(90.0, color=colors["neutral"], linestyle="--", linewidth=0.9, label="Nominal 90%")
    ax.set(xlabel="Oscillator period (s)", ylabel="Empirical coverage (%)", title="a  Spatial-block interval coverage")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.set_ylim(75, 100)
    ax.grid(alpha=0.2, linewidth=0.5)
    ax.legend(frameon=False, loc="lower left")

    ax = axes[0, 1]
    fold_metrics = coverage[coverage["scope"].eq("fold")]
    grouped = fold_metrics.groupby("period_s")["median_interval_factor_span"]
    span_median = grouped.median().reindex(overall["period_s"])
    span_min = grouped.min().reindex(overall["period_s"])
    span_max = grouped.max().reindex(overall["period_s"])
    ax.fill_between(
        overall["period_s"],
        span_min,
        span_max,
        color=colors["interval"],
        alpha=0.22,
        linewidth=0,
        label="Range across five held blocks",
    )
    ax.semilogx(
        overall["period_s"],
        span_median,
        marker="o",
        color=colors["interval"],
        label="Median block",
    )
    ax.set(
        xlabel="Oscillator period (s)",
        ylabel="Upper/lower multiplier factor",
        title="b  Width of the calibrated interval",
    )
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.grid(alpha=0.2, linewidth=0.5)
    ax.legend(frameon=False, loc="upper right")

    sa3 = intervals[np.isclose(intervals["period_s"], 3.0)].copy()
    ax = axes[1, 0]
    sorted_sa3 = sa3.sort_values("point_prediction_log10").reset_index(drop=True)
    selected_indices = np.linspace(0, len(sorted_sa3) - 1, 35).round().astype(int)
    selected = sorted_sa3.iloc[selected_indices].copy()
    x = np.arange(len(selected))
    ax.vlines(
        x,
        selected["lower_multiplier"],
        selected["upper_multiplier"],
        color="#A0A0A0",
        linewidth=1.1,
    )
    ax.scatter(x, selected["point_multiplier"], color=colors["point"], s=13, zorder=3)
    ax.axhline(1.0, color=colors["neutral"], linestyle="--", linewidth=0.8)
    ax.set_yscale("log")
    ax.set(
        xlabel="Stations ordered by point multiplier",
        ylabel="SA(3.0 s) station multiplier",
        title="c  Prediction intervals across station terms",
    )
    ax.set_xticks([])
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)

    ax = axes[1, 1]
    city = city_cases[city_cases["city"].eq("Tokyo")].copy()
    if city.empty:
        city = city_cases[city_cases["city"].eq(city_cases["city"].iloc[0])].copy()
    siteid2 = int(city["siteid2"].iloc[0])
    probability = "50y_10pct"
    spectrum = hazard_intervals[
        hazard_intervals["siteid2"].eq(siteid2)
        & hazard_intervals["probability_level"].eq(probability)
    ].sort_values("period_s")
    city_name = str(city["city"].iloc[0])
    ax.fill_between(
        spectrum["period_s"],
        spectrum["adjusted_surface_sa_lower_g"],
        spectrum["adjusted_surface_sa_upper_g"],
        color=colors["interval"],
        alpha=0.22,
        linewidth=0,
        label="Station-model 90% interval",
    )
    ax.loglog(
        spectrum["period_s"],
        spectrum["ergodic_surface_sa_g"],
        color=colors["neutral"],
        marker="o",
        linewidth=1.2,
        label="Ergodic surface reference",
    )
    ax.loglog(
        spectrum["period_s"],
        spectrum["adjusted_surface_sa_point_g"],
        color=colors["point"],
        marker="s",
        linewidth=1.3,
        label="Station-adjusted point estimate",
    )
    ax.set(
        xlabel="Oscillator period (s)",
        ylabel="Spectral acceleration (g)",
        title=f"d  {city_name} matched site, 50 yr 10%",
    )
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.grid(alpha=0.2, which="both", linewidth=0.5)
    ax.legend(frameon=False, loc="best")

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    intervals: pd.DataFrame,
    coverage: pd.DataFrame,
    hazard_intervals: pd.DataFrame,
    lower_quantile: float,
    upper_quantile: float,
) -> None:
    overall = coverage[coverage["scope"].eq("overall")].sort_values("period_s")
    sa3 = overall[np.isclose(overall["period_s"], 3.0)].iloc[0]
    sa3_intervals = intervals[np.isclose(intervals["period_s"], 3.0)]
    nominal = 100.0 * (upper_quantile - lower_quantile)
    lines = [
        "# Spatial-block station-term uncertainty propagation",
        "",
        "## Design",
        "",
        f"- Target interval: weighted empirical {nominal:.0f}% prediction interval ({lower_quantile:.2f}--{upper_quantile:.2f} error quantiles).",
        "- Calibration: for each held spatial fold, signed out-of-fold errors are estimated only from the other folds.",
        "- Weighting: station errors are weighted by the number of contributing strong-motion records.",
        "- Propagation: lower, point, and upper station terms multiply the same matched J-SHIS surface-reference spectrum.",
        "",
        "## Results",
        "",
        f"- Station-period intervals: {len(intervals):,} rows.",
        f"- Propagated station-period-probability values: {len(hazard_intervals):,} rows.",
        f"- SA(3.0 s) station coverage: {sa3['station_coverage_pct']:.1f}%.",
        f"- SA(3.0 s) record-weighted coverage: {sa3['record_weighted_coverage_pct']:.1f}%.",
        f"- SA(3.0 s) median upper/lower multiplier span: {sa3_intervals['interval_factor_span'].median():.2f}.",
        f"- Coverage range across eight periods: {overall['station_coverage_pct'].min():.1f}%--{overall['station_coverage_pct'].max():.1f}%.",
        "",
        "## Interpretation boundary",
        "",
        "These intervals describe prediction error for the station-term model under the existing spatial-block validation design. They do not include uncertainty in earthquake occurrence, source characterization, the MF2013 backbone, J-SHIS hazard calculations, or the surface-reference conversion. The propagated bounds are conditional sensitivity coordinates, not complete non-ergodic PSHA confidence intervals.",
        "",
    ]
    OUT_AUDIT.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lower-quantile", type=float, default=0.05)
    parser.add_argument("--upper-quantile", type=float, default=0.95)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.lower_quantile < args.upper_quantile < 1.0:
        raise ValueError("quantiles must satisfy 0 < lower < upper < 1")
    predictions = pd.read_csv(PREDICTIONS, low_memory=False)
    hazard_values = pd.read_csv(HAZARD_VALUES, low_memory=False)
    city_cases = pd.read_csv(CITY_CASES, low_memory=False)
    intervals, coverage = calibrate_intervals(predictions, args.lower_quantile, args.upper_quantile)
    hazard_intervals, hazard_summary = propagate_to_spectra(hazard_values, intervals)

    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    intervals.to_csv(OUT_INTERVALS, index=False)
    coverage.to_csv(OUT_COVERAGE, index=False)
    hazard_intervals.to_csv(OUT_HAZARD_INTERVALS, index=False)
    hazard_summary.to_csv(OUT_HAZARD_SUMMARY, index=False)
    plot_figure(intervals, coverage, hazard_intervals, city_cases)
    write_audit(intervals, coverage, hazard_intervals, args.lower_quantile, args.upper_quantile)
    print(f"wrote {OUT_AUDIT}", flush=True)


if __name__ == "__main__":
    main()
