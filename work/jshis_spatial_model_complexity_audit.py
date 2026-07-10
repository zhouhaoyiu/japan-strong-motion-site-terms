#!/usr/bin/env python3
"""Test spatial station prediction across fixed model-complexity settings."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core


OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_spatial_model_complexity_sensitivity.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_spatial_model_complexity_sensitivity.md"

SETTINGS: dict[str, dict[str, float | int]] = {
    "primary": {},
    "fewer_leaves": {"model__max_leaf_nodes": 7},
    "more_leaves": {"model__max_leaf_nodes": 31},
    "larger_minimum_leaf": {"model__min_samples_leaf": 70},
    "stronger_l2": {"model__l2_regularization": 1.0},
    "shorter_iteration_budget": {"model__max_iter": 200},
}


def evaluate(args: argparse.Namespace) -> pd.DataFrame:
    station_terms = pd.read_csv(args.station_terms)
    selected = station_terms[
        station_terms["model"].eq("mf2013_site")
        & station_terms["n_records"].ge(args.min_records)
        & station_terms["lon"].notna()
        & station_terms["lat"].notna()
    ].copy()
    rows: list[dict[str, float | int | str]] = []
    for spec in core.PERIODS:
        frame = selected[selected["period_s"].eq(spec.period_s)].copy().reset_index(drop=True)
        target = frame["station_effect_log10"].to_numpy(float)
        weights = frame["n_records"].to_numpy(float)
        folds = core.make_spatial_folds(frame, n_splits=args.spatial_folds, seed=args.seed)
        for setting, overrides in SETTINGS.items():
            out_of_fold = np.full(len(frame), np.nan, dtype=float)
            for fold in sorted(np.unique(folds)):
                train = folds != fold
                test = folds == fold
                estimator = core.build_model(include_space=True, seed=args.seed + int(fold))
                if overrides:
                    estimator.set_params(**overrides)
                estimator.fit(
                    frame.loc[train],
                    target[train],
                    model__sample_weight=weights[train],
                )
                out_of_fold[test] = estimator.predict(frame.loc[test])
            centering_shift = float(np.average(out_of_fold, weights=weights))
            prediction = out_of_fold - centering_shift
            for fold in sorted(np.unique(folds)):
                test = folds == fold
                baseline_rmse = core.weighted_rmse(
                    target[test], np.zeros(test.sum()), weights[test]
                )
                rmse = core.weighted_rmse(target[test], prediction[test], weights[test])
                rows.append(
                    {
                        "period_s": spec.period_s,
                        "period_code": spec.code,
                        "setting": setting,
                        "spatial_seed": args.seed,
                        "scope": "fold",
                        "fold": int(fold),
                        "n_stations": int(test.sum()),
                        "n_records_weight": int(weights[test].sum()),
                        "weighted_rmse_log10": rmse,
                        "zero_baseline_rmse_log10": baseline_rmse,
                        "rmse_reduction_vs_zero_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse,
                        "observed_predicted_correlation": float(
                            np.corrcoef(target[test], prediction[test])[0, 1]
                        ),
                        "prediction_centering_shift_log10": centering_shift,
                    }
                )
            baseline_rmse = core.weighted_rmse(target, np.zeros_like(target), weights)
            rmse = core.weighted_rmse(target, prediction, weights)
            rows.append(
                {
                    "period_s": spec.period_s,
                    "period_code": spec.code,
                    "setting": setting,
                    "spatial_seed": args.seed,
                    "scope": "overall",
                    "fold": -1,
                    "n_stations": len(frame),
                    "n_records_weight": int(weights.sum()),
                    "weighted_rmse_log10": rmse,
                    "zero_baseline_rmse_log10": baseline_rmse,
                    "rmse_reduction_vs_zero_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse,
                    "observed_predicted_correlation": float(np.corrcoef(target, prediction)[0, 1]),
                    "prediction_centering_shift_log10": centering_shift,
                }
            )
        print(f"complexity audit: period={spec.period_s:g}s", flush=True)
    return pd.DataFrame(rows)


def write_audit(metrics: pd.DataFrame) -> None:
    overall = metrics[metrics["scope"].eq("overall")]
    sa3 = overall[overall["period_s"].eq(3.0)].sort_values("setting")
    primary = overall[overall["setting"].eq("primary")]
    lines = [
        "# Spatial station-model complexity sensitivity",
        "",
        f"The six settings use the same station targets, spatial blocks, preprocessing, weights and random seed "
        f"({int(metrics['spatial_seed'].iloc[0])}). "
        "Only tree complexity, regularization or iteration budget changes. The primary setting remains fixed; "
        "the alternatives are a sensitivity audit and are not used to select a replacement model.",
        "",
        "| Setting | SA(3.0 s) RMSE gain (%) | SA(3.0 s) correlation |",
        "|---|---:|---:|",
    ]
    for row in sa3.itertuples(index=False):
        lines.append(
            f"| {row.setting} | {row.rmse_reduction_vs_zero_pct:.2f} | "
            f"{row.observed_predicted_correlation:.3f} |"
        )
    lines.extend(
        [
            "",
            f"Across all periods and settings, overall RMSE gains range from "
            f"{overall['rmse_reduction_vs_zero_pct'].min():.2f}% to "
            f"{overall['rmse_reduction_vs_zero_pct'].max():.2f}%.",
            f"The primary setting reproduces eight period rows and its minimum gain is "
            f"{primary['rmse_reduction_vs_zero_pct'].min():.2f}%.",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    core.SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    metrics = evaluate(args)
    metrics.to_csv(OUT_METRICS, index=False)
    write_audit(metrics)
    print(f"wrote {OUT_METRICS}")
    print(f"wrote {OUT_AUDIT}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--station-terms", type=Path, default=core.OUT_STATION_TERMS)
    parser.add_argument("--min-records", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
