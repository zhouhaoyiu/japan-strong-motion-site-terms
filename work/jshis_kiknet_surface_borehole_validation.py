#!/usr/bin/env python3
"""Validate surface station terms with paired KiK-net borehole spectra."""

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
from jshis_robustness_stress_tests import weighted_correlation


OUT_STATIONS = core.SUPPLEMENT_DIR / "jshis_kiknet_surface_borehole_station_transfer.csv"
OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_kiknet_surface_borehole_validation.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_kiknet_surface_borehole_validation.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "figure_kiknet_surface_borehole_validation.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "figure_kiknet_surface_borehole_validation.png"


def load_pairs(flatfile: Path, site: pd.DataFrame) -> pd.DataFrame:
    response_columns = [spec.rotd100_col for spec in core.PERIODS]
    with zipfile.ZipFile(flatfile) as archive:
        records = pd.read_csv(
            archive.open("smrec_schema.tsv"),
            sep="\t",
            usecols=["smrec_id", "site_id", "eq_source_id", *response_columns],
            low_memory=False,
        )
    records["siteid2"] = pd.to_numeric(records["site_id"], errors="coerce").astype("Int64") // 10
    metadata = site[
        [
            "siteid2",
            "site_code",
            "lon",
            "lat",
            "sensor_depth_glminus",
            "obs_network_id",
            "installation_situation_id",
        ]
    ]
    records = records.merge(metadata, on="siteid2", how="left", validate="many_to_one")
    records = records[
        records["obs_network_id"].eq(2)
        & records["installation_situation_id"].isin([1, 2])
    ].copy()

    keys = ["eq_source_id", "site_code"]
    surface = records[records["installation_situation_id"].eq(1)].copy()
    borehole = records[records["installation_situation_id"].eq(2)].copy()
    if surface.duplicated(keys).any() or borehole.duplicated(keys).any():
        raise ValueError("KiK-net event-site records are not unique")

    surface = surface.rename(
        columns={
            "siteid2": "surface_siteid2",
            "sensor_depth_glminus": "surface_sensor_depth_m",
            **{column: f"{column}_surface" for column in response_columns},
        }
    )
    borehole = borehole.rename(
        columns={
            "siteid2": "borehole_siteid2",
            "sensor_depth_glminus": "borehole_sensor_depth_m",
            **{column: f"{column}_borehole" for column in response_columns},
        }
    )
    return surface[
        keys
        + [
            "surface_siteid2",
            "lon",
            "lat",
            "surface_sensor_depth_m",
            *[f"{column}_surface" for column in response_columns],
        ]
    ].merge(
        borehole[
            keys
            + [
                "borehole_siteid2",
                "borehole_sensor_depth_m",
                *[f"{column}_borehole" for column in response_columns],
            ]
        ],
        on=keys,
        how="inner",
        validate="one_to_one",
    )


def period_pairs(pairs: pd.DataFrame, spec: core.PeriodSpec) -> pd.DataFrame:
    surface = pd.to_numeric(pairs[f"{spec.rotd100_col}_surface"], errors="coerce")
    borehole = pd.to_numeric(pairs[f"{spec.rotd100_col}_borehole"], errors="coerce")
    usable = surface.gt(0) & borehole.gt(0) & np.isfinite(surface) & np.isfinite(borehole)
    frame = pairs.loc[
        usable,
        [
            "eq_source_id",
            "site_code",
            "surface_siteid2",
            "borehole_siteid2",
            "lon",
            "lat",
            "borehole_sensor_depth_m",
        ],
    ].copy()
    frame["log10_surface_borehole_ratio"] = np.log10(
        surface.loc[usable].to_numpy() / borehole.loc[usable].to_numpy()
    )
    return frame


def station_summary(frame: pd.DataFrame, spec: core.PeriodSpec) -> pd.DataFrame:
    summary = (
        frame.groupby(
            [
                "surface_siteid2",
                "borehole_siteid2",
                "site_code",
                "lon",
                "lat",
                "borehole_sensor_depth_m",
            ],
            as_index=False,
        )["log10_surface_borehole_ratio"]
        .agg(
            n_paired_records="size",
            transfer_log10_median="median",
            transfer_log10_q25=lambda values: values.quantile(0.25),
            transfer_log10_q75=lambda values: values.quantile(0.75),
        )
    )
    summary["period_s"] = spec.period_s
    summary["period_code"] = spec.code
    return summary


def correlations(frame: pd.DataFrame) -> tuple[float, float, float]:
    x = frame["station_effect_log10"].to_numpy(float)
    y = frame["transfer_log10_median"].to_numpy(float)
    weights = np.minimum(
        frame["n_records"].to_numpy(float), frame["n_paired_records"].to_numpy(float)
    )
    return (
        float(np.corrcoef(x, y)[0, 1]),
        float(pd.Series(x).corr(pd.Series(y), method="spearman")),
        weighted_correlation(x, y, weights),
    )


def depth_adjusted_correlation(frame: pd.DataFrame) -> float:
    """Pearson correlation after removing linear log-borehole-depth trends."""
    depth = pd.to_numeric(frame["borehole_sensor_depth_m"], errors="coerce")
    usable = depth.gt(0) & frame["station_effect_log10"].notna() & frame[
        "transfer_log10_median"
    ].notna()
    clean = frame.loc[usable]
    log_depth = np.log10(depth.loc[usable].to_numpy(float))
    design = np.column_stack([np.ones(len(clean)), log_depth])
    station = clean["station_effect_log10"].to_numpy(float)
    transfer = clean["transfer_log10_median"].to_numpy(float)
    station_residual = station - design @ np.linalg.lstsq(design, station, rcond=None)[0]
    transfer_residual = transfer - design @ np.linalg.lstsq(design, transfer, rcond=None)[0]
    return float(np.corrcoef(station_residual, transfer_residual)[0, 1])


def fit_weighted_line(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    design = np.column_stack([np.ones(len(x)), x])
    root_weight = np.sqrt(weights)[:, None]
    intercept, slope = np.linalg.lstsq(design * root_weight, y * root_weight[:, 0], rcond=None)[0]
    return float(intercept), float(slope)


def event_folds(events: np.ndarray, n_splits: int, seed: int) -> dict[int, int]:
    shuffled = np.random.default_rng(seed).permutation(np.sort(events))
    return {int(event_id): index % n_splits for index, event_id in enumerate(shuffled)}


def validate_period(
    residuals: pd.DataFrame,
    paired: pd.DataFrame,
    spec: core.PeriodSpec,
    n_splits: int,
    seed: int,
    min_train_records: int,
    min_test_records: int,
) -> list[dict[str, float | int | str]]:
    fold_map = event_folds(residuals["eq_source_id"].dropna().unique(), n_splits, seed)
    residuals = residuals.copy()
    residuals["fold"] = residuals["eq_source_id"].map(fold_map).astype(int)
    paired = paired.copy()
    paired["fold"] = paired["eq_source_id"].map(fold_map)
    paired = paired.dropna(subset="fold")
    paired["fold"] = paired["fold"].astype(int)

    rows: list[dict[str, float | int | str]] = []
    for fold in range(n_splits):
        train_records = residuals[residuals["fold"].ne(fold)]
        test_records = residuals[residuals["fold"].eq(fold)]
        train_terms, _, train_fit = core.two_way_decomposition(train_records, "residual_site")
        test_terms, _, test_fit = core.two_way_decomposition(test_records, "residual_site")
        train_terms = train_terms.rename(
            columns={"station_effect_log10": "train_station_term", "n_records": "n_train_records"}
        )
        test_terms = test_terms.rename(
            columns={"station_effect_log10": "test_station_term", "n_records": "n_test_records"}
        )

        train_transfer = (
            paired[paired["fold"].ne(fold)]
            .groupby("surface_siteid2")["log10_surface_borehole_ratio"]
            .agg(train_transfer="median", n_train_pairs="size")
            .reset_index()
        )
        test_transfer = (
            paired[paired["fold"].eq(fold)]
            .groupby("surface_siteid2")["log10_surface_borehole_ratio"]
            .agg(test_transfer="median", n_test_pairs="size")
            .reset_index()
        )
        repeated = train_transfer.merge(test_transfer, on="surface_siteid2", validate="one_to_one")
        repeated = repeated[
            repeated["n_train_pairs"].ge(min_train_records)
            & repeated["n_test_pairs"].ge(min_test_records)
        ]

        calibration = train_terms.merge(
            train_transfer, left_on="siteid2", right_on="surface_siteid2", validate="one_to_one"
        )
        calibration = calibration[
            calibration["n_train_records"].ge(min_train_records)
            & calibration["n_train_pairs"].ge(min_train_records)
        ]
        calibration_weights = np.minimum(
            calibration["n_train_records"].to_numpy(float),
            calibration["n_train_pairs"].to_numpy(float),
        )
        intercept, slope = fit_weighted_line(
            calibration["train_transfer"].to_numpy(float),
            calibration["train_station_term"].to_numpy(float),
            calibration_weights,
        )

        validation = test_terms.merge(
            train_transfer, left_on="siteid2", right_on="surface_siteid2", validate="one_to_one"
        )
        validation = validation[
            validation["n_test_records"].ge(min_test_records)
            & validation["n_train_pairs"].ge(min_train_records)
        ]
        weights = validation["n_test_records"].to_numpy(float)
        observed = np.array(validation["test_station_term"].to_numpy(float), copy=True)
        predicted = np.array(
            intercept + slope * validation["train_transfer"].to_numpy(float), copy=True
        )
        observed -= np.average(observed, weights=weights)
        predicted -= np.average(predicted, weights=weights)
        baseline_rmse = core.weighted_rmse(observed, np.zeros_like(observed), weights)
        prediction_rmse = core.weighted_rmse(observed, predicted, weights)
        rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "scope": "fold",
                "fold": fold,
                "event_split_seed": seed,
                "n_train_events": train_records["eq_source_id"].nunique(),
                "n_test_events": test_records["eq_source_id"].nunique(),
                "n_transfer_repeat_stations": len(repeated),
                "transfer_train_test_pearson": repeated["train_transfer"].corr(
                    repeated["test_transfer"]
                ),
                "transfer_train_test_spearman": repeated["train_transfer"].corr(
                    repeated["test_transfer"], method="spearman"
                ),
                "n_calibration_stations": len(calibration),
                "n_validation_stations": len(validation),
                "calibration_intercept": intercept,
                "calibration_slope": slope,
                "test_station_term_pearson": float(np.corrcoef(observed, predicted)[0, 1]),
                "test_station_term_weighted_pearson": weighted_correlation(
                    observed, predicted, weights
                ),
                "zero_baseline_rmse_log10": baseline_rmse,
                "transfer_prediction_rmse_log10": prediction_rmse,
                "transfer_prediction_rmse_gain_pct": 100.0
                * (baseline_rmse - prediction_rmse)
                / baseline_rmse,
                "train_solver_max_change": train_fit["max_parameter_change"],
                "test_solver_max_change": test_fit["max_parameter_change"],
            }
        )

    fold_frame = pd.DataFrame(rows)
    mean_row: dict[str, float | int | str] = {
        "period_s": spec.period_s,
        "period_code": spec.code,
        "scope": "fold_mean",
        "fold": -1,
        "event_split_seed": seed,
    }
    for column in fold_frame.columns:
        if column not in mean_row and column not in {"period_code", "scope"}:
            mean_row[column] = float(fold_frame[column].mean())
    rows.append(mean_row)
    return rows


def plot_results(stations: pd.DataFrame, metrics: pd.DataFrame) -> None:
    full = metrics[metrics["scope"].eq("full")].sort_values("period_s")
    means = metrics[metrics["scope"].eq("fold_mean")].sort_values("period_s")
    folds = metrics[metrics["scope"].eq("fold")]
    sa3 = stations[stations["period_s"].eq(3.0)].copy()
    linked = sa3.dropna(subset="station_effect_log10")

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), constrained_layout=True)
    ax = axes[0, 0]
    points = ax.scatter(
        sa3["lon"],
        sa3["lat"],
        c=sa3["transfer_log10_median"],
        cmap="RdBu_r",
        vmin=-0.5,
        vmax=0.5,
        s=13,
        linewidths=0,
    )
    fig.colorbar(points, ax=ax, label=r"Median $\log_{10}$(surface/borehole)")
    ax.set(xlabel="Longitude (°E)", ylabel="Latitude (°N)", title="a  KiK-net paired transfer at 3.0 s")

    ax = axes[0, 1]
    ax.scatter(
        linked["transfer_log10_median"],
        linked["station_effect_log10"],
        s=np.clip(linked["n_paired_records"], 8, 80),
        alpha=0.55,
        color="#2878B5",
        linewidths=0,
    )
    slope, intercept = np.polyfit(
        linked["transfer_log10_median"], linked["station_effect_log10"], 1
    )
    xline = np.linspace(linked["transfer_log10_median"].min(), linked["transfer_log10_median"].max(), 100)
    ax.plot(xline, intercept + slope * xline, color="#C43C39", lw=1.6)
    row3 = full[full["period_s"].eq(3.0)].iloc[0]
    ax.text(
        0.04,
        0.96,
        f"Pearson $r$ = {row3['station_transfer_pearson']:.3f}\nSpearman $\\rho$ = {row3['station_transfer_spearman']:.3f}",
        transform=ax.transAxes,
        va="top",
    )
    ax.set(
        xlabel=r"Median $\log_{10}$(surface/borehole)",
        ylabel="MF2013 station term (log$_{10}$)",
        title="b  Paired-sensor correspondence",
    )

    ax = axes[1, 0]
    ax.plot(full["period_s"], full["station_transfer_pearson"], "o-", label="Station field, Pearson")
    ax.plot(full["period_s"], full["station_transfer_spearman"], "s-", label="Station field, Spearman")
    ax.plot(means["period_s"], means["transfer_train_test_pearson"], "^-", label="Transfer across events")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.5", "1", "2", "3", "5"])
    ax.set(xlabel="Period (s)", ylabel="Correlation", ylim=(0, 1.02), title="c  Period dependence and event stability")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    periods = means["period_s"].to_numpy(float)
    correlations = means["test_station_term_pearson"].to_numpy(float)
    corr_min = folds.groupby("period_s")["test_station_term_pearson"].min().reindex(periods).to_numpy()
    corr_max = folds.groupby("period_s")["test_station_term_pearson"].max().reindex(periods).to_numpy()
    ax.fill_between(periods, corr_min, corr_max, color="#2878B5", alpha=0.15)
    ax.plot(periods, correlations, "o-", color="#2878B5", label="Correlation")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.5", "1", "2", "3", "5"])
    ax.set(xlabel="Period (s)", ylabel="Test-event correlation", ylim=(0, 0.8), title="d  Held-event prediction from paired transfer")
    twin = ax.twinx()
    twin.plot(
        periods,
        means["transfer_prediction_rmse_gain_pct"],
        "s--",
        color="#C43C39",
        label="RMSE gain",
    )
    twin.set_ylabel("RMSE gain (%)", color="#C43C39")
    twin.tick_params(axis="y", colors="#C43C39")
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = twin.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, frameon=False, fontsize=8, loc="upper right")

    for axis in axes.flat:
        axis.grid(alpha=0.18, linewidth=0.6)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    coefficients = core.load_coefficients(args.coefficients)
    site, source = core.load_metadata(args.flatfile)
    pairs = load_pairs(args.flatfile, site)
    residuals = core.load_record_residuals(
        args.flatfile, coefficients, site, source, chunksize=args.chunksize
    )
    primary_terms = pd.read_csv(args.station_terms)
    primary_terms = primary_terms[
        primary_terms["model"].eq("mf2013_site")
        & primary_terms["n_records"].ge(args.minimum_station_records)
    ]

    station_blocks = []
    metric_rows: list[dict[str, float | int | str]] = []
    for spec in core.PERIODS:
        paired = period_pairs(pairs, spec)
        stations = station_summary(paired, spec)
        terms = primary_terms[primary_terms["period_s"].eq(spec.period_s)][
            ["siteid2", "n_records", "station_effect_log10"]
        ]
        stations = stations.merge(
            terms, left_on="surface_siteid2", right_on="siteid2", how="left", validate="one_to_one"
        ).drop(columns="siteid2")
        linked = stations.dropna(subset="station_effect_log10")
        pearson, spearman, weighted = correlations(linked)
        depth_adjusted = depth_adjusted_correlation(linked)
        station_blocks.append(stations)
        metric_rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "scope": "full",
                "fold": -1,
                "event_split_seed": args.seed,
                "n_paired_records": len(paired),
                "n_paired_events": paired["eq_source_id"].nunique(),
                "n_paired_stations": stations["surface_siteid2"].nunique(),
                "n_transfer_stations_min_records": stations[
                    "n_paired_records"
                ].ge(args.minimum_station_records).sum(),
                "n_linked_station_terms": len(linked),
                "station_transfer_pearson": pearson,
                "station_transfer_spearman": spearman,
                "station_transfer_weighted_pearson": weighted,
                "station_transfer_depth_adjusted_pearson": depth_adjusted,
            }
        )
        metric_rows.extend(
            validate_period(
                residuals[spec.period_s],
                paired,
                spec,
                n_splits=args.folds,
                seed=args.seed,
                min_train_records=args.minimum_train_records,
                min_test_records=args.minimum_test_records,
            )
        )
        print(f"KiK-net validation: period={spec.period_s:g}s", flush=True)

    station_frame = pd.concat(station_blocks, ignore_index=True)
    metric_frame = pd.DataFrame(metric_rows)
    station_frame.to_csv(OUT_STATIONS, index=False)
    metric_frame.to_csv(OUT_METRICS, index=False)
    plot_results(station_frame, metric_frame)

    full3 = metric_frame[metric_frame["scope"].eq("full") & metric_frame["period_s"].eq(3.0)].iloc[0]
    mean3 = metric_frame[
        metric_frame["scope"].eq("fold_mean") & metric_frame["period_s"].eq(3.0)
    ].iloc[0]
    lines = [
        "# KiK-net surface/borehole spectral validation",
        "",
        "Surface and borehole RotD100 spectra are paired by earthquake and KiK-net site code. Transfer ratios are estimated from disjoint training events and used to predict station terms fitted to held-out surface records.",
        "",
        f"- Paired records at each period: {int(full3['n_paired_records']):,}",
        f"- Paired sites and earthquakes: {int(full3['n_paired_stations']):,} and {int(full3['n_paired_events']):,}",
        f"- SA(3.0 s) station-term/transfer Pearson and Spearman correlations: {full3['station_transfer_pearson']:.3f} and {full3['station_transfer_spearman']:.3f}",
        f"- SA(3.0 s) Pearson correlation after log-borehole-depth adjustment: {full3['station_transfer_depth_adjusted_pearson']:.3f}",
        f"- SA(3.0 s) transfer train/test correlation across event folds: {mean3['transfer_train_test_pearson']:.3f}",
        f"- SA(3.0 s) cross-event station-term correlation: {mean3['test_station_term_pearson']:.3f}",
        f"- SA(3.0 s) cross-event RMSE gain over a zero station term: {mean3['transfer_prediction_rmse_gain_pct']:.1f}%",
        "",
        "The paired-sensor result provides an independent physical correlate of the event-adjusted surface station field. It is a validation of station-level response, not a calibrated nonlinear transfer-function model.",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT_STATIONS}", flush=True)
    print(f"wrote {OUT_METRICS}", flush=True)
    print(f"wrote {OUT_AUDIT}", flush=True)
    print(f"wrote {OUT_FIGURE_PDF}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=core.DEFAULT_COEFFICIENTS)
    parser.add_argument("--station-terms", type=Path, default=core.OUT_STATION_TERMS)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20_260_710)
    parser.add_argument("--minimum-station-records", type=int, default=20)
    parser.add_argument("--minimum-train-records", type=int, default=15)
    parser.add_argument("--minimum-test-records", type=int, default=5)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
