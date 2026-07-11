#!/usr/bin/env python3
"""Validate frozen station-term models across networks and unseen events."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

import jshis_event_adjusted_station_model as core
from jshis_robustness_stress_tests import weighted_correlation


OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_cross_network_transfer_metrics.csv"
OUT_PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_cross_network_transfer_predictions.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_cross_network_transfer_validation.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "figure_cross_network_transfer_validation.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "figure_cross_network_transfer_validation.png"


def make_model(include_space: bool, seed: int) -> Pipeline:
    features = core.PHYSICAL_FEATURES + (core.SPATIAL_FEATURES if include_space else [])
    preprocess = ColumnTransformer(
        [
            (
                "numeric",
                SimpleImputer(
                    strategy="median",
                    add_indicator=True,
                    keep_empty_features=True,
                ),
                features,
            )
        ]
    )
    estimator = HistGradientBoostingRegressor(
        loss="squared_error",
        learning_rate=0.04,
        max_iter=400,
        max_leaf_nodes=15,
        min_samples_leaf=35,
        l2_regularization=0.1,
        random_state=seed,
    )
    return Pipeline([("preprocess", preprocess), ("model", estimator)])


def assign_event_folds(events: np.ndarray, n_splits: int, seed: int) -> dict[int, int]:
    shuffled = np.random.default_rng(seed).permutation(np.sort(events.astype(int)))
    return {int(event): int(index % n_splits) for index, event in enumerate(shuffled)}


def supported_common_events(frame: pd.DataFrame, minimum_records: int) -> np.ndarray:
    counts = (
        frame.groupby(["eq_source_id", "network_label"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["K-NET", "KiK-net"], fill_value=0)
    )
    return counts.index[
        counts["K-NET"].ge(minimum_records)
        & counts["KiK-net"].ge(minimum_records)
    ].to_numpy(int)


def station_terms(
    frame: pd.DataFrame,
    site: pd.DataFrame,
    network: str,
    minimum_records: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    selected = frame[frame["network_label"].eq(network)]
    terms, _, fit = core.two_way_decomposition(selected, "residual_site")
    terms = terms[terms["n_records"].ge(minimum_records)].copy()
    metadata = site.drop_duplicates("siteid2")
    terms = terms.merge(metadata, on="siteid2", how="left", validate="one_to_one")
    return terms, fit


def metric_values(frame: pd.DataFrame) -> dict[str, float]:
    observed = frame["observed_station_term_log10"].to_numpy(float)
    predicted = frame["frozen_prediction_log10"].to_numpy(float)
    weights = frame["n_test_records"].to_numpy(float)
    baseline_rmse = core.weighted_rmse(observed, np.zeros_like(observed), weights)
    prediction_rmse = core.weighted_rmse(observed, predicted, weights)
    baseline_mae = core.weighted_mae(observed, np.zeros_like(observed), weights)
    prediction_mae = core.weighted_mae(observed, predicted, weights)
    return {
        "pearson": float(np.corrcoef(observed, predicted)[0, 1]),
        "weighted_pearson": weighted_correlation(observed, predicted, weights),
        "zero_baseline_rmse_log10": baseline_rmse,
        "prediction_rmse_log10": prediction_rmse,
        "rmse_gain_pct": 100.0 * (baseline_rmse - prediction_rmse) / baseline_rmse,
        "zero_baseline_mae_log10": baseline_mae,
        "prediction_mae_log10": prediction_mae,
        "mae_gain_pct": 100.0 * (baseline_mae - prediction_mae) / baseline_mae,
        "prediction_weighted_mean_log10": float(np.average(predicted, weights=weights)),
    }


def bootstrap_intervals(
    frame: pd.DataFrame,
    replicates: int,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    values: dict[str, list[float]] = {
        "pearson": [],
        "weighted_pearson": [],
        "rmse_gain_pct": [],
        "mae_gain_pct": [],
    }
    size = len(frame)
    for _ in range(replicates):
        sample = frame.iloc[rng.integers(0, size, size=size)]
        metrics = metric_values(sample)
        for key in values:
            values[key].append(metrics[key])
    intervals: dict[str, float] = {}
    for key, samples in values.items():
        lower, upper = np.quantile(samples, [0.025, 0.975])
        intervals[f"{key}_ci_low"] = float(lower)
        intervals[f"{key}_ci_high"] = float(upper)
    return intervals


def aggregate_stations(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    identity = [
        "period_s",
        "period_code",
        "source_network",
        "target_network",
        "model",
        "siteid2",
        "site_code",
        "lon",
        "lat",
        "avs30",
        "d1400",
    ]
    for keys, block in predictions.groupby(identity, dropna=False, sort=False):
        weights = block["n_test_records"].to_numpy(float)
        row = dict(zip(identity, keys, strict=True))
        row.update(
            {
                "n_test_records": int(weights.sum()),
                "n_test_folds": int(block["fold"].nunique()),
                "observed_station_term_log10": float(
                    np.average(block["observed_station_term_log10"], weights=weights)
                ),
                "frozen_prediction_log10": float(
                    np.average(block["frozen_prediction_log10"], weights=weights)
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def run_direction(
    frame: pd.DataFrame,
    site: pd.DataFrame,
    spec: core.PeriodSpec,
    fold_map: dict[int, int],
    source_network: str,
    target_network: str,
    args: argparse.Namespace,
) -> tuple[list[dict[str, float | int | str]], list[pd.DataFrame]]:
    frame = frame.copy()
    frame["fold"] = frame["eq_source_id"].map(fold_map)
    frame = frame.dropna(subset="fold")
    frame["fold"] = frame["fold"].astype(int)
    metrics: list[dict[str, float | int | str]] = []
    predictions: list[pd.DataFrame] = []
    for fold in range(args.folds):
        train_records = frame[frame["fold"].ne(fold)]
        test_records = frame[frame["fold"].eq(fold)]
        source_terms, source_fit = station_terms(
            train_records,
            site,
            source_network,
            args.minimum_train_records,
        )
        target_terms, target_fit = station_terms(
            test_records,
            site,
            target_network,
            args.minimum_test_records,
        )
        for model_name, include_space in [
            ("physical_hgb", False),
            ("physical_spatial_hgb", True),
        ]:
            model = make_model(include_space, args.seed + fold)
            train_weights = source_terms["n_records"].to_numpy(float)
            model.fit(
                source_terms,
                source_terms["station_effect_log10"].to_numpy(float),
                model__sample_weight=train_weights,
            )
            train_prediction = model.predict(source_terms)
            centering_shift = float(np.average(train_prediction, weights=train_weights))
            predicted = model.predict(target_terms) - centering_shift
            block = target_terms[
                ["siteid2", "site_code", "lon", "lat", "avs30", "d1400", "n_records"]
            ].copy()
            block = block.rename(
                columns={
                    "n_records": "n_test_records",
                }
            )
            block["period_s"] = spec.period_s
            block["period_code"] = spec.code
            block["source_network"] = source_network
            block["target_network"] = target_network
            block["model"] = model_name
            block["fold"] = fold
            block["observed_station_term_log10"] = target_terms[
                "station_effect_log10"
            ].to_numpy(float)
            block["frozen_prediction_log10"] = predicted
            block["prediction_error_log10"] = (
                block["observed_station_term_log10"] - block["frozen_prediction_log10"]
            )
            predictions.append(block)

            row: dict[str, float | int | str] = {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "source_network": source_network,
                "target_network": target_network,
                "model": model_name,
                "scope": "fold",
                "fold": fold,
                "event_split_seed": args.seed,
                "n_train_events": int(
                    train_records.loc[
                        train_records["network_label"].eq(source_network), "eq_source_id"
                    ].nunique()
                ),
                "n_test_events": int(
                    test_records.loc[
                        test_records["network_label"].eq(target_network), "eq_source_id"
                    ].nunique()
                ),
                "n_train_stations": len(source_terms),
                "n_test_stations": len(target_terms),
                "n_train_records": int(source_terms["n_records"].sum()),
                "n_test_records": int(target_terms["n_records"].sum()),
                "training_prediction_center_log10": centering_shift,
                "train_solver_max_change": source_fit["max_parameter_change"],
                "test_solver_max_change": target_fit["max_parameter_change"],
            }
            row.update(metric_values(block))
            metrics.append(row)
    return metrics, predictions


def summarize_aggregate(
    aggregate: pd.DataFrame,
    fold_metrics: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    rows = []
    group_columns = [
        "period_s",
        "period_code",
        "source_network",
        "target_network",
        "model",
    ]
    for keys, block in aggregate.groupby(group_columns, sort=False):
        period_s, period_code, source_network, target_network, model_name = keys
        matching_folds = fold_metrics[
            fold_metrics["period_s"].eq(period_s)
            & fold_metrics["source_network"].eq(source_network)
            & fold_metrics["target_network"].eq(target_network)
            & fold_metrics["model"].eq(model_name)
        ]
        row: dict[str, float | int | str] = {
            "period_s": period_s,
            "period_code": period_code,
            "source_network": source_network,
            "target_network": target_network,
            "model": model_name,
            "scope": "station_aggregate",
            "fold": -1,
            "event_split_seed": args.seed,
            "n_train_events": int(round(matching_folds["n_train_events"].mean())),
            "n_test_events": int(round(matching_folds["n_test_events"].mean())),
            "n_train_stations": int(round(matching_folds["n_train_stations"].mean())),
            "n_test_stations": len(block),
            "n_train_records": int(round(matching_folds["n_train_records"].mean())),
            "n_test_records": int(block["n_test_records"].sum()),
            "training_prediction_center_log10": float(
                matching_folds["training_prediction_center_log10"].mean()
            ),
            "train_solver_max_change": float(
                matching_folds["train_solver_max_change"].max()
            ),
            "test_solver_max_change": float(
                matching_folds["test_solver_max_change"].max()
            ),
        }
        row.update(metric_values(block))
        row.update(
            bootstrap_intervals(
                block,
                replicates=args.bootstrap_replicates,
                seed=args.seed + int(round(float(period_s) * 1000)) + len(rows),
            )
        )
        rows.append(row)
    return pd.DataFrame(rows)


def plot_results(aggregate: pd.DataFrame, metrics: pd.DataFrame) -> None:
    selected = metrics[
        metrics["scope"].eq("station_aggregate")
        & metrics["source_network"].eq("K-NET")
        & metrics["target_network"].eq("KiK-net")
    ].sort_values(["model", "period_s"])
    prediction3 = aggregate[
        aggregate["period_s"].eq(3.0)
        & aggregate["source_network"].eq("K-NET")
        & aggregate["target_network"].eq("KiK-net")
        & aggregate["model"].eq("physical_spatial_hgb")
    ].copy()
    row3 = selected[
        selected["period_s"].eq(3.0)
        & selected["model"].eq("physical_spatial_hgb")
    ].iloc[0]

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), constrained_layout=True)
    ax = axes[0, 0]
    points = ax.scatter(
        prediction3["lon"],
        prediction3["lat"],
        c=prediction3["observed_station_term_log10"],
        cmap="RdBu_r",
        vmin=-0.45,
        vmax=0.45,
        s=16,
        linewidths=0,
    )
    fig.colorbar(points, ax=ax, label="Held-event station term (log$_{10}$)")
    ax.set(
        xlabel="Longitude (°E)",
        ylabel="Latitude (°N)",
        title="a  Unseen KiK-net stations and events at 3.0 s",
    )

    ax = axes[0, 1]
    weights = prediction3["n_test_records"].to_numpy(float)
    ax.scatter(
        prediction3["frozen_prediction_log10"],
        prediction3["observed_station_term_log10"],
        s=np.clip(weights, 8, 75),
        alpha=0.5,
        color="#2878B5",
        linewidths=0,
    )
    limits = np.quantile(
        np.r_[
            prediction3["frozen_prediction_log10"],
            prediction3["observed_station_term_log10"],
        ],
        [0.005, 0.995],
    )
    ax.plot(limits, limits, color="#444444", lw=1.2, ls="--")
    ax.text(
        0.04,
        0.96,
        (
            f"Pearson $r$ = {row3['pearson']:.3f} "
            f"[{row3['pearson_ci_low']:.3f}, {row3['pearson_ci_high']:.3f}]\n"
            f"RMSE gain = {row3['rmse_gain_pct']:.1f}%"
        ),
        transform=ax.transAxes,
        va="top",
    )
    ax.set(
        xlim=limits,
        ylim=limits,
        xlabel="Frozen K-NET prediction (log$_{10}$)",
        ylabel="Observed KiK-net station term (log$_{10}$)",
        title="b  Network- and event-external prediction",
    )

    colors = {"physical_hgb": "#D55E00", "physical_spatial_hgb": "#0072B2"}
    labels = {"physical_hgb": "Physical", "physical_spatial_hgb": "Physical + spatial"}
    ax = axes[1, 0]
    for model_name, block in selected.groupby("model"):
        block = block.sort_values("period_s")
        ax.fill_between(
            block["period_s"],
            block["pearson_ci_low"],
            block["pearson_ci_high"],
            color=colors[model_name],
            alpha=0.13,
        )
        ax.plot(
            block["period_s"],
            block["pearson"],
            "o-",
            color=colors[model_name],
            label=labels[model_name],
        )
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.5", "1", "2", "3", "5"])
    ax.set(
        xlabel="Period (s)",
        ylabel="Observed-predicted correlation",
        title="c  Transfer across response periods",
    )
    ax.legend(frameon=False)

    ax = axes[1, 1]
    for model_name, block in selected.groupby("model"):
        block = block.sort_values("period_s")
        ax.fill_between(
            block["period_s"],
            block["rmse_gain_pct_ci_low"],
            block["rmse_gain_pct_ci_high"],
            color=colors[model_name],
            alpha=0.13,
        )
        ax.plot(
            block["period_s"],
            block["rmse_gain_pct"],
            "s-",
            color=colors[model_name],
            label=labels[model_name],
        )
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.5", "1", "2", "3", "5"])
    ax.set(
        xlabel="Period (s)",
        ylabel="RMSE gain over zero term (%)",
        title="d  Out-of-network error reduction",
    )
    ax.legend(frameon=False)

    for axis in axes.flat:
        axis.grid(alpha=0.18, linewidth=0.6)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    coefficients = core.load_coefficients(args.coefficients)
    site, source = core.load_metadata(args.flatfile)
    records = core.load_record_residuals(
        args.flatfile,
        coefficients,
        site,
        source,
        chunksize=args.chunksize,
    )
    site_network = site[["siteid2", "network_label"]].drop_duplicates("siteid2")

    metric_rows: list[dict[str, float | int | str]] = []
    prediction_blocks: list[pd.DataFrame] = []
    for spec in core.PERIODS:
        frame = records[spec.period_s].merge(
            site_network,
            on="siteid2",
            how="left",
            validate="many_to_one",
        )
        common_events = supported_common_events(frame, args.minimum_event_network_records)
        fold_map = assign_event_folds(common_events, args.folds, args.seed)
        frame = frame[frame["eq_source_id"].isin(common_events)].copy()
        for source_network, target_network in [
            ("K-NET", "KiK-net"),
            ("KiK-net", "K-NET"),
        ]:
            direction_metrics, direction_predictions = run_direction(
                frame,
                site,
                spec,
                fold_map,
                source_network,
                target_network,
                args,
            )
            metric_rows.extend(direction_metrics)
            prediction_blocks.extend(direction_predictions)
        print(
            f"cross-network: period={spec.period_s:g}s common_events={len(common_events):,}",
            flush=True,
        )

    predictions = pd.concat(prediction_blocks, ignore_index=True)
    aggregate = aggregate_stations(predictions)
    fold_metrics = pd.DataFrame(metric_rows)
    aggregate_metrics = summarize_aggregate(aggregate, fold_metrics, args)
    metrics = pd.concat([fold_metrics, aggregate_metrics], ignore_index=True, sort=False)
    predictions.to_csv(OUT_PREDICTIONS, index=False)
    metrics.to_csv(OUT_METRICS, index=False)
    plot_results(aggregate, metrics)

    primary = metrics[
        metrics["scope"].eq("station_aggregate")
        & metrics["source_network"].eq("K-NET")
        & metrics["target_network"].eq("KiK-net")
        & metrics["model"].eq("physical_spatial_hgb")
        & metrics["period_s"].eq(3.0)
    ].iloc[0]
    reverse = metrics[
        metrics["scope"].eq("station_aggregate")
        & metrics["source_network"].eq("KiK-net")
        & metrics["target_network"].eq("K-NET")
        & metrics["model"].eq("physical_spatial_hgb")
        & metrics["period_s"].eq(3.0)
    ].iloc[0]
    lines = [
        "# Cross-network and held-event station-term validation",
        "",
        "K-NET station terms are estimated only from training events. The fixed public-variable model is then evaluated against KiK-net station terms estimated only from disjoint held-out events. The reverse direction is retained as a transfer stress test.",
        "",
        f"- Event folds: {args.folds}; fixed seed: {args.seed}.",
        f"- SA(3.0 s) K-NET to KiK-net target stations: {int(primary['n_test_stations']):,}.",
        f"- SA(3.0 s) Pearson correlation: {primary['pearson']:.3f} (station bootstrap 95% CI {primary['pearson_ci_low']:.3f} to {primary['pearson_ci_high']:.3f}).",
        f"- SA(3.0 s) RMSE gain: {primary['rmse_gain_pct']:.1f}% (95% CI {primary['rmse_gain_pct_ci_low']:.1f}% to {primary['rmse_gain_pct_ci_high']:.1f}%).",
        f"- Reverse KiK-net to K-NET correlation and RMSE gain: {reverse['pearson']:.3f} and {reverse['rmse_gain_pct']:.1f}%.",
        "",
        "The validation separates both recording network and earthquake set. No target-network station term, target event, or target prediction is used to fit or recenter the frozen model.",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT_METRICS}", flush=True)
    print(f"wrote {OUT_PREDICTIONS}", flush=True)
    print(f"wrote {OUT_AUDIT}", flush=True)
    print(f"wrote {OUT_FIGURE_PDF}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=core.DEFAULT_COEFFICIENTS)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20_260_711)
    parser.add_argument("--minimum-event-network-records", type=int, default=5)
    parser.add_argument("--minimum-train-records", type=int, default=15)
    parser.add_argument("--minimum-test-records", type=int, default=5)
    parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
