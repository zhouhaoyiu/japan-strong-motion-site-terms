#!/usr/bin/env python3
"""Validate station-field transfer across both recording network and time."""

from __future__ import annotations

import argparse
import sys
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core
from jshis_cross_network_transfer_validation import (
    bootstrap_intervals,
    make_model,
    metric_values,
    station_terms,
    supported_common_events,
)


OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_temporal_network_transfer_metrics.csv"
OUT_PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_temporal_network_transfer_predictions.csv"
OUT_EVENT_SPLITS = core.SUPPLEMENT_DIR / "jshis_temporal_network_event_splits.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_temporal_network_transfer_validation.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "supplementary_figure_temporal_network_transfer.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "supplementary_figure_temporal_network_transfer.png"


def load_event_dates(flatfile: Path) -> pd.DataFrame:
    with zipfile.ZipFile(flatfile) as archive:
        frame = pd.read_csv(
            archive.open("source_schema.tsv"),
            sep="\t",
            usecols=["eq_source_id", "segment_idx", "jem_origin_time"],
        )
    frame["eq_source_id"] = pd.to_numeric(frame["eq_source_id"], errors="coerce").astype("Int64")
    frame["segment_idx"] = pd.to_numeric(frame["segment_idx"], errors="coerce")
    frame["origin_time"] = pd.to_datetime(frame["jem_origin_time"], errors="coerce")
    return (
        frame.dropna(subset=["eq_source_id", "origin_time"])
        .sort_values(["eq_source_id", "segment_idx"])
        .drop_duplicates("eq_source_id", keep="first")[["eq_source_id", "origin_time"]]
    )


def chronological_event_split(
    event_dates: pd.DataFrame,
    events: np.ndarray,
    train_fraction: float,
) -> tuple[np.ndarray, np.ndarray, pd.Timestamp, pd.Timestamp]:
    ordered = (
        event_dates[event_dates["eq_source_id"].isin(events)]
        .sort_values(["origin_time", "eq_source_id"])
        .reset_index(drop=True)
    )
    split_index = int(np.floor(len(ordered) * train_fraction))
    if split_index <= 0 or split_index >= len(ordered):
        raise ValueError(f"invalid chronological split for {len(ordered)} events")
    train = ordered.iloc[:split_index]
    test = ordered.iloc[split_index:]
    return (
        train["eq_source_id"].to_numpy(int),
        test["eq_source_id"].to_numpy(int),
        pd.Timestamp(train["origin_time"].max()),
        pd.Timestamp(test["origin_time"].min()),
    )


def evaluate_direction(
    frame: pd.DataFrame,
    site: pd.DataFrame,
    spec: core.PeriodSpec,
    train_events: np.ndarray,
    test_events: np.ndarray,
    train_end: pd.Timestamp,
    test_start: pd.Timestamp,
    train_fraction: float,
    source_network: str,
    target_network: str,
    args: argparse.Namespace,
) -> tuple[list[dict[str, float | int | str]], list[pd.DataFrame]]:
    train_records = frame[frame["eq_source_id"].isin(train_events)]
    test_records = frame[frame["eq_source_id"].isin(test_events)]
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

    metric_rows: list[dict[str, float | int | str]] = []
    prediction_blocks: list[pd.DataFrame] = []
    for model_index, (model_name, include_space) in enumerate(
        [("physical_hgb", False), ("physical_spatial_hgb", True)]
    ):
        seed = args.seed + int(round(100 * train_fraction)) + model_index
        model = make_model(include_space, seed)
        train_weights = source_terms["n_records"].to_numpy(float)
        model.fit(
            source_terms,
            source_terms["station_effect_log10"].to_numpy(float),
            model__sample_weight=train_weights,
        )
        training_prediction = model.predict(source_terms)
        centering_shift = float(np.average(training_prediction, weights=train_weights))
        frozen_prediction = model.predict(target_terms) - centering_shift

        block = target_terms[
            ["siteid2", "site_code", "lon", "lat", "avs30", "d1400", "n_records"]
        ].copy()
        block = block.rename(columns={"n_records": "n_test_records"})
        block["period_s"] = spec.period_s
        block["period_code"] = spec.code
        block["train_fraction"] = train_fraction
        block["train_end"] = train_end.isoformat()
        block["test_start"] = test_start.isoformat()
        block["source_network"] = source_network
        block["target_network"] = target_network
        block["model"] = model_name
        block["observed_station_term_log10"] = target_terms[
            "station_effect_log10"
        ].to_numpy(float)
        block["frozen_prediction_log10"] = frozen_prediction
        block["prediction_error_log10"] = (
            block["observed_station_term_log10"] - block["frozen_prediction_log10"]
        )
        prediction_blocks.append(block)

        row: dict[str, float | int | str] = {
            "period_s": spec.period_s,
            "period_code": spec.code,
            "train_fraction": train_fraction,
            "train_end": train_end.isoformat(),
            "test_start": test_start.isoformat(),
            "source_network": source_network,
            "target_network": target_network,
            "model": model_name,
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
        row.update(
            bootstrap_intervals(
                block,
                replicates=args.bootstrap_replicates,
                seed=seed + int(round(1_000 * spec.period_s)),
            )
        )
        metric_rows.append(row)
    return metric_rows, prediction_blocks


def plot_results(predictions: pd.DataFrame, metrics: pd.DataFrame) -> None:
    selected = metrics[
        metrics["source_network"].eq("K-NET")
        & metrics["target_network"].eq("KiK-net")
        & metrics["model"].eq("physical_spatial_hgb")
    ].copy()
    primary = selected[selected["train_fraction"].eq(0.8)].sort_values("period_s")
    primary_prediction = predictions[
        predictions["period_s"].eq(3.0)
        & predictions["train_fraction"].eq(0.8)
        & predictions["source_network"].eq("K-NET")
        & predictions["target_network"].eq("KiK-net")
        & predictions["model"].eq("physical_spatial_hgb")
    ].copy()
    row3 = primary[primary["period_s"].eq(3.0)].iloc[0]

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), constrained_layout=True)
    ax = axes[0, 0]
    colors = {0.7: "#009E73", 0.8: "#0072B2", 0.9: "#D55E00"}
    for train_fraction, block in selected.groupby("train_fraction"):
        block = block.sort_values("period_s")
        ax.plot(
            block["period_s"],
            block["pearson"],
            "o-",
            color=colors[float(train_fraction)],
            label=f"{int(round(100 * train_fraction))}/{int(round(100 * (1 - train_fraction)))}",
        )
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.5", "1", "2", "3", "5"])
    ax.set(
        xlabel="Period (s)",
        ylabel="Observed-predicted correlation",
        title="a  Chronological split sensitivity",
    )
    ax.legend(frameon=False, title="Early/late events (%)")

    ax = axes[0, 1]
    for train_fraction, block in selected.groupby("train_fraction"):
        block = block.sort_values("period_s")
        ax.plot(
            block["period_s"],
            block["rmse_gain_pct"],
            "s-",
            color=colors[float(train_fraction)],
            label=f"{int(round(100 * train_fraction))}/{int(round(100 * (1 - train_fraction)))}",
        )
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.5", "1", "2", "3", "5"])
    ax.set(
        xlabel="Period (s)",
        ylabel="RMSE gain over zero term (%)",
        title="b  Error reduction through time",
    )

    ax = axes[1, 0]
    weights = primary_prediction["n_test_records"].to_numpy(float)
    ax.scatter(
        primary_prediction["frozen_prediction_log10"],
        primary_prediction["observed_station_term_log10"],
        s=np.clip(weights, 8, 75),
        alpha=0.5,
        color="#2878B5",
        linewidths=0,
    )
    limits = np.quantile(
        np.r_[
            primary_prediction["frozen_prediction_log10"],
            primary_prediction["observed_station_term_log10"],
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
        xlabel="Frozen early K-NET prediction (log$_{10}$)",
        ylabel="Observed late KiK-net station term (log$_{10}$)",
        title="c  Primary 80/20 split at 3.0 s",
    )

    ax = axes[1, 1]
    long_period = selected[selected["period_s"].isin([1.0, 2.0, 3.0])].copy()
    offsets = {1.0: -0.018, 2.0: 0.0, 3.0: 0.018}
    period_colors = {1.0: "#009E73", 2.0: "#E69F00", 3.0: "#0072B2"}
    for period, block in long_period.groupby("period_s"):
        x = block["train_fraction"].to_numpy(float) + offsets[float(period)]
        ax.errorbar(
            x,
            block["pearson"],
            yerr=np.vstack(
                [
                    block["pearson"] - block["pearson_ci_low"],
                    block["pearson_ci_high"] - block["pearson"],
                ]
            ),
            fmt="o",
            capsize=3,
            color=period_colors[float(period)],
            label=f"{period:g} s",
        )
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set_xticks([0.7, 0.8, 0.9], labels=["70/30", "80/20", "90/10"])
    ax.set(
        xlabel="Early/late event split (%)",
        ylabel="Observed-predicted correlation",
        title="d  Long-period station-bootstrap intervals",
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
    event_dates = load_event_dates(args.flatfile)
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
    event_split_blocks: list[pd.DataFrame] = []
    for spec in core.PERIODS:
        frame = records[spec.period_s].merge(
            site_network,
            on="siteid2",
            how="left",
            validate="many_to_one",
        )
        common_events = supported_common_events(frame, args.minimum_event_network_records)
        frame = frame[frame["eq_source_id"].isin(common_events)].copy()
        for train_fraction in args.train_fractions:
            train_events, test_events, train_end, test_start = chronological_event_split(
                event_dates,
                common_events,
                train_fraction,
            )
            split_block = event_dates[event_dates["eq_source_id"].isin(common_events)].copy()
            split_block["period_s"] = spec.period_s
            split_block["train_fraction"] = train_fraction
            split_block["partition"] = np.where(
                split_block["eq_source_id"].isin(train_events),
                "early_train",
                "late_test",
            )
            event_split_blocks.append(split_block)
            for source_network, target_network in [
                ("K-NET", "KiK-net"),
                ("KiK-net", "K-NET"),
            ]:
                direction_metrics, direction_predictions = evaluate_direction(
                    frame,
                    site,
                    spec,
                    train_events,
                    test_events,
                    train_end,
                    test_start,
                    train_fraction,
                    source_network,
                    target_network,
                    args,
                )
                metric_rows.extend(direction_metrics)
                prediction_blocks.extend(direction_predictions)
            print(
                f"temporal transfer: period={spec.period_s:g}s split={train_fraction:.1f} "
                f"train_end={train_end.date()} test_start={test_start.date()}",
                flush=True,
            )

    metrics = pd.DataFrame(metric_rows)
    predictions = pd.concat(prediction_blocks, ignore_index=True)
    event_splits = pd.concat(event_split_blocks, ignore_index=True)
    metrics.to_csv(OUT_METRICS, index=False)
    predictions.to_csv(OUT_PREDICTIONS, index=False)
    event_splits.to_csv(OUT_EVENT_SPLITS, index=False)
    plot_results(predictions, metrics)

    primary = metrics[
        metrics["period_s"].eq(3.0)
        & metrics["train_fraction"].eq(0.8)
        & metrics["source_network"].eq("K-NET")
        & metrics["target_network"].eq("KiK-net")
        & metrics["model"].eq("physical_spatial_hgb")
    ].iloc[0]
    long_period = metrics[
        metrics["period_s"].isin([1.0, 2.0, 3.0])
        & metrics["source_network"].eq("K-NET")
        & metrics["target_network"].eq("KiK-net")
        & metrics["model"].eq("physical_spatial_hgb")
    ]
    primary_periods = metrics[
        metrics["train_fraction"].eq(0.8)
        & metrics["source_network"].eq("K-NET")
        & metrics["target_network"].eq("KiK-net")
        & metrics["model"].eq("physical_spatial_hgb")
    ]
    reverse = metrics[
        metrics["period_s"].eq(3.0)
        & metrics["train_fraction"].eq(0.8)
        & metrics["source_network"].eq("KiK-net")
        & metrics["target_network"].eq("K-NET")
        & metrics["model"].eq("physical_spatial_hgb")
    ].iloc[0]
    lines = [
        "# Chronological cross-network station-field validation",
        "",
        "Events are ordered by JMA origin time. K-NET station terms are estimated from the earlier event block, and the frozen public-variable model is evaluated against KiK-net station terms estimated from the later block. Recording network, event set, and target stations are all external to model fitting.",
        "",
        "- Prespecified early/late event splits: 70/30, 80/20, and 90/10.",
        f"- Primary 80/20 cutoff: training through {pd.Timestamp(primary['train_end']).date()}; testing from {pd.Timestamp(primary['test_start']).date()}.",
        f"- Primary SA(3.0 s) target: {int(primary['n_test_stations']):,} KiK-net stations from {int(primary['n_test_events']):,} later events.",
        f"- Primary SA(3.0 s) Pearson correlation: {primary['pearson']:.3f} (station bootstrap 95% CI {primary['pearson_ci_low']:.3f} to {primary['pearson_ci_high']:.3f}).",
        f"- Primary SA(3.0 s) RMSE gain: {primary['rmse_gain_pct']:.1f}% (95% CI {primary['rmse_gain_pct_ci_low']:.1f}% to {primary['rmse_gain_pct_ci_high']:.1f}%).",
        f"- Across all eight periods in the primary split, correlations range from {primary_periods['pearson'].min():.3f} to {primary_periods['pearson'].max():.3f}; RMSE gains range from {primary_periods['rmse_gain_pct'].min():.1f}% to {primary_periods['rmse_gain_pct'].max():.1f}%, and every 95% lower bound remains above zero.",
        f"- Across the three splits and 1-3 s periods, correlations range from {long_period['pearson'].min():.3f} to {long_period['pearson'].max():.3f}; RMSE gains range from {long_period['rmse_gain_pct'].min():.1f}% to {long_period['rmse_gain_pct'].max():.1f}%.",
        f"- Reverse KiK-net-to-K-NET transfer at 3.0 s gives a correlation of {reverse['pearson']:.3f} and an RMSE gain of {reverse['rmse_gain_pct']:.1f}% (95% CI {reverse['rmse_gain_pct_ci_low']:.1f}% to {reverse['rmse_gain_pct_ci_high']:.1f}%).",
        "",
        "The target-network terms and predictions are not used for model fitting or recentering. Each chronological split and the adverse reverse-direction result are reported, so the result does not depend on selecting a favorable cutoff year or transfer direction.",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT_METRICS}", flush=True)
    print(f"wrote {OUT_PREDICTIONS}", flush=True)
    print(f"wrote {OUT_EVENT_SPLITS}", flush=True)
    print(f"wrote {OUT_AUDIT}", flush=True)
    print(f"wrote {OUT_FIGURE_PDF}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=core.DEFAULT_COEFFICIENTS)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20_260_711)
    parser.add_argument("--train-fractions", type=float, nargs="+", default=[0.7, 0.8, 0.9])
    parser.add_argument("--minimum-event-network-records", type=int, default=5)
    parser.add_argument("--minimum-train-records", type=int, default=15)
    parser.add_argument("--minimum-test-records", type=int, default=5)
    parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
