#!/usr/bin/env python3
"""Audit station-field results under selected published MF2013 regression screens."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core
from jshis_robustness_stress_tests import weighted_correlation, weighted_quantile


OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_mf2013_applicability_sensitivity.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_mf2013_applicability_sensitivity.md"


def restricted_records(
    records: dict[float, pd.DataFrame],
    source: pd.DataFrame,
    minimum_mw: float,
    maximum_distance_km: float,
    minimum_event_stations: int,
) -> dict[float, pd.DataFrame]:
    source_magnitude = source[["eq_source_id", "mw"]]
    selected: dict[float, pd.DataFrame] = {}
    for period_s, frame in records.items():
        joined = frame.merge(source_magnitude, on="eq_source_id", how="left", validate="many_to_one")
        subset = joined[
            joined["mw"].ge(minimum_mw) & joined["fault_dist"].lt(maximum_distance_km)
        ].copy()
        event_station_counts = subset.groupby("eq_source_id")["siteid2"].nunique()
        retained_events = event_station_counts[event_station_counts.ge(minimum_event_stations)].index
        selected[period_s] = subset[subset["eq_source_id"].isin(retained_events)].drop(columns="mw")
    return selected


def build_station_terms(
    records: dict[float, pd.DataFrame],
    site: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    term_blocks = []
    fit_rows = []
    for spec in core.PERIODS:
        frame = records[spec.period_s]
        terms, _, metrics = core.two_way_decomposition(frame, "residual_site")
        terms = terms.merge(site, on="siteid2", how="left", validate="one_to_one")
        terms["period_s"] = spec.period_s
        terms["period_code"] = spec.code
        terms["model"] = "mf2013_site"
        terms["response_component"] = "RotD100"
        term_blocks.append(terms)
        fit_rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "n_records": len(frame),
                "n_events": frame["eq_source_id"].nunique(),
                "n_stations": frame["siteid2"].nunique(),
                "solver_iterations": int(metrics["iterations"]),
                "solver_max_parameter_change": float(metrics["max_parameter_change"]),
            }
        )
    return pd.concat(term_blocks, ignore_index=True), pd.DataFrame(fit_rows)


def compare_station_fields(
    restricted_terms: pd.DataFrame,
    full_terms: pd.DataFrame,
    minimum_records: int,
) -> pd.DataFrame:
    rows = []
    for spec in core.PERIODS:
        restricted = restricted_terms[
            restricted_terms["period_s"].eq(spec.period_s)
            & restricted_terms["n_records"].ge(minimum_records)
        ][["siteid2", "n_records", "station_effect_log10"]]
        full = full_terms[
            full_terms["period_s"].eq(spec.period_s)
            & full_terms["model"].eq("mf2013_site")
            & full_terms["n_records"].ge(minimum_records)
        ][["siteid2", "n_records", "station_effect_log10"]]
        paired = full.merge(restricted, on="siteid2", suffixes=("_full", "_restricted"))
        weights = paired[["n_records_full", "n_records_restricted"]].min(axis=1).to_numpy(float)
        full_values = paired["station_effect_log10_full"].to_numpy(float, copy=True)
        restricted_values = paired["station_effect_log10_restricted"].to_numpy(float, copy=True)
        full_values -= np.average(full_values, weights=weights)
        restricted_values -= np.average(restricted_values, weights=weights)
        difference = restricted_values - full_values
        rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "n_eligible_stations": len(restricted),
                "n_paired_stations": len(paired),
                "station_field_correlation": weighted_correlation(
                    full_values, restricted_values, weights
                ),
                "station_field_difference_q95_log10": weighted_quantile(
                    np.abs(difference), 0.95, weights
                ),
            }
        )
    return pd.DataFrame(rows)


def compare_predictions(
    restricted_predictions: pd.DataFrame,
    full_predictions: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for spec in core.PERIODS:
        restricted = restricted_predictions[
            restricted_predictions["period_s"].eq(spec.period_s)
            & restricted_predictions["model"].eq("physical_spatial_hgb")
        ][["siteid2", "n_records", "centered_oof_prediction_log10"]]
        full = full_predictions[
            full_predictions["period_s"].eq(spec.period_s)
            & full_predictions["residual_model"].eq("mf2013_site")
            & full_predictions["model"].eq("physical_spatial_hgb")
        ][["siteid2", "n_records", "centered_oof_prediction_log10"]]
        paired = full.merge(restricted, on="siteid2", suffixes=("_full", "_restricted"))
        weights = paired[["n_records_full", "n_records_restricted"]].min(axis=1).to_numpy(float)
        full_values = paired["centered_oof_prediction_log10_full"].to_numpy(float, copy=True)
        restricted_values = paired["centered_oof_prediction_log10_restricted"].to_numpy(
            float, copy=True
        )
        full_values -= np.average(full_values, weights=weights)
        restricted_values -= np.average(restricted_values, weights=weights)
        rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "prediction_paired_stations": len(paired),
                "oof_prediction_correlation": weighted_correlation(
                    full_values, restricted_values, weights
                ),
                "oof_prediction_difference_q95_log10": weighted_quantile(
                    np.abs(restricted_values - full_values), 0.95, weights
                ),
            }
        )
    return pd.DataFrame(rows)


def compare_hazard_values(
    restricted_values: pd.DataFrame,
    full_values: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for spec in core.PERIODS:
        restricted = restricted_values[
            restricted_values["period_s"].eq(spec.period_s)
            & restricted_values["probability_level"].eq("50y_10pct")
        ][
            [
                "siteid2",
                "ergodic_surface_sa_g",
                "oof_multiplier",
                "adjusted_surface_sa_g",
            ]
        ]
        full = full_values[
            full_values["period_s"].eq(spec.period_s)
            & full_values["probability_level"].eq("50y_10pct")
        ][["siteid2", "ergodic_surface_sa_g", "oof_multiplier", "adjusted_surface_sa_g"]]
        paired = full.merge(restricted, on="siteid2", suffixes=("_primary", "_restricted"))
        rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "common_hazard_stations": len(paired),
                "common_surface_sa_median_g": paired[
                    "ergodic_surface_sa_g_restricted"
                ].median(),
                "common_primary_multiplier_q05": paired["oof_multiplier_primary"].quantile(0.05),
                "common_primary_multiplier_q50": paired["oof_multiplier_primary"].quantile(0.50),
                "common_primary_multiplier_q95": paired["oof_multiplier_primary"].quantile(0.95),
                "common_primary_adjusted_sa_median_g": paired[
                    "adjusted_surface_sa_g_primary"
                ].median(),
                "common_restricted_adjusted_sa_median_g": paired[
                    "adjusted_surface_sa_g_restricted"
                ].median(),
            }
        )
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> None:
    coefficients = core.load_coefficients(args.coefficients)
    site, source = core.load_metadata(args.flatfile)
    records = core.load_record_residuals(
        args.flatfile, coefficients, site, source, chunksize=args.chunksize
    )
    restricted = restricted_records(
        records,
        source,
        minimum_mw=args.minimum_mw,
        maximum_distance_km=args.maximum_distance_km,
        minimum_event_stations=args.minimum_event_stations,
    )
    station_terms, fits = build_station_terms(restricted, site)
    full_terms = pd.read_csv(args.full_station_terms)
    field_comparison = compare_station_fields(station_terms, full_terms, args.minimum_station_records)

    event_holdout = core.event_holdout_repeatability(
        restricted, n_splits=args.spatial_folds, seed=args.seed
    )
    event_summary = event_holdout[event_holdout["scope"].eq("fold_mean")][
        ["period_s", "train_test_station_correlation", "rmse_reduction_vs_zero_pct"]
    ].rename(
        columns={
            "train_test_station_correlation": "restricted_event_holdout_correlation",
            "rmse_reduction_vs_zero_pct": "restricted_event_holdout_rmse_gain_pct",
        }
    )

    predictions, model_metrics = core.cross_validate_station_models(
        station_terms,
        min_records=args.minimum_station_records,
        n_splits=args.spatial_folds,
        seed=args.seed,
    )
    spatial_summary = model_metrics[
        model_metrics["model"].eq("physical_spatial_hgb")
        & model_metrics["scope"].eq("overall")
    ][["period_s", "rmse_reduction_vs_zero_pct", "observed_predicted_correlation"]].rename(
        columns={
            "rmse_reduction_vs_zero_pct": "restricted_spatial_rmse_gain_pct",
            "observed_predicted_correlation": "restricted_spatial_correlation",
        }
    )
    full_predictions = pd.read_csv(args.full_predictions)
    prediction_comparison = compare_predictions(predictions, full_predictions)

    hazard_values, hazard_summary = core.build_hazard_values(
        predictions, coefficients, args.response_map
    )
    hazard_10pct = hazard_summary[hazard_summary["probability_level"].eq("50y_10pct")][
        [
            "period_s",
            "n_station_values",
            "official_vs400_sa_g_q50",
            "ergodic_surface_sa_g_q50",
            "oof_multiplier_q05",
            "oof_multiplier_q50",
            "oof_multiplier_q95",
            "adjusted_surface_sa_g_q50",
        ]
    ].rename(
        columns={
            "n_station_values": "restricted_hazard_stations",
            "official_vs400_sa_g_q50": "restricted_official_sa_median_g",
            "ergodic_surface_sa_g_q50": "restricted_surface_sa_median_g",
            "oof_multiplier_q05": "restricted_multiplier_q05",
            "oof_multiplier_q50": "restricted_multiplier_q50",
            "oof_multiplier_q95": "restricted_multiplier_q95",
            "adjusted_surface_sa_g_q50": "restricted_adjusted_sa_median_g",
        }
    )
    full_hazard_values = pd.read_csv(args.full_hazard_values)
    common_hazard = compare_hazard_values(hazard_values, full_hazard_values)

    metrics = (
        fits.merge(field_comparison, on=["period_s", "period_code"], validate="one_to_one")
        .merge(event_summary, on="period_s", validate="one_to_one")
        .merge(spatial_summary, on="period_s", validate="one_to_one")
        .merge(prediction_comparison, on=["period_s", "period_code"], validate="one_to_one")
        .merge(hazard_10pct, on="period_s", validate="one_to_one")
        .merge(common_hazard, on=["period_s", "period_code"], validate="one_to_one")
    )
    metrics.insert(2, "minimum_mw", args.minimum_mw)
    metrics.insert(3, "maximum_distance_km_exclusive", args.maximum_distance_km)
    metrics.insert(4, "minimum_event_stations", args.minimum_event_stations)
    metrics.to_csv(OUT_METRICS, index=False)

    sa3 = metrics[metrics["period_s"].eq(3.0)].iloc[0]
    lines = [
        "# MF2013 magnitude, distance and event-station sensitivity",
        "",
        "Morikawa and Fujiwara (2013) selected regression records with Mw >= 5.5, ground-surface sensors, at least five triggered stations per event and source distance below 200 km. This sensitivity applies the magnitude, distance and event-station-count conditions to the primary surface sample. It does not reconstruct the paper's additional magnitude-dependent PGA truncation.",
        "",
        f"- Records at each period: {int(sa3['n_records']):,}",
        f"- Earthquakes: {int(sa3['n_events']):,}",
        f"- Stations with at least {args.minimum_station_records} records: {int(sa3['n_eligible_stations']):,}",
        f"- SA(3.0 s) station-field correlation with the primary field: {sa3['station_field_correlation']:.3f}",
        f"- SA(3.0 s) station-field 95th-percentile absolute difference: {sa3['station_field_difference_q95_log10']:.3f} log10 units",
        f"- SA(3.0 s) restricted event-holdout correlation: {sa3['restricted_event_holdout_correlation']:.3f}",
        f"- SA(3.0 s) restricted spatial RMSE gain: {sa3['restricted_spatial_rmse_gain_pct']:.1f}%",
        f"- SA(3.0 s) primary/restricted OOF-prediction correlation: {sa3['oof_prediction_correlation']:.3f}",
        f"- SA(3.0 s) restricted multiplier 5th--95th percentiles: {sa3['restricted_multiplier_q05']:.3f}--{sa3['restricted_multiplier_q95']:.3f}",
        f"- SA(3.0 s) restricted surface and adjusted medians: {sa3['restricted_surface_sa_median_g']:.3f} and {sa3['restricted_adjusted_sa_median_g']:.3f} g",
        f"- On the same {int(sa3['common_hazard_stations']):,} stations, the primary adjusted median is {sa3['common_primary_adjusted_sa_median_g']:.3f} g and the restricted adjusted median is {sa3['common_restricted_adjusted_sa_median_g']:.3f} g, from a common surface median of {sa3['common_surface_sa_median_g']:.3f} g",
        "",
        "The restricted calculation tests selected magnitude, distance and event-station conditions from the original regression design. The primary analysis retains the broader public-flatfile domain because it evaluates the implementation used with the national response-spectrum product. Results outside the implemented screens are interpreted through this sensitivity.",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False), flush=True)
    print(f"wrote {OUT_METRICS}", flush=True)
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=core.DEFAULT_COEFFICIENTS)
    parser.add_argument("--response-map", type=Path, default=core.DEFAULT_RESPONSE_MAP)
    parser.add_argument("--full-station-terms", type=Path, default=core.OUT_STATION_TERMS)
    parser.add_argument("--full-predictions", type=Path, default=core.OUT_MODEL_PREDICTIONS)
    parser.add_argument("--full-hazard-values", type=Path, default=core.OUT_HAZARD_VALUES)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--minimum-mw", type=float, default=5.5)
    parser.add_argument("--maximum-distance-km", type=float, default=200.0)
    parser.add_argument("--minimum-event-stations", type=int, default=5)
    parser.add_argument("--minimum-station-records", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
