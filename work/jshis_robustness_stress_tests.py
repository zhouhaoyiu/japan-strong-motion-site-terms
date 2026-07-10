#!/usr/bin/env python3
"""Stress-test station-term conclusions across networks, regions, and events."""

from __future__ import annotations

import argparse
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


SUPPLEMENT_DIR = core.SUPPLEMENT_DIR
FIGURE_DIR = core.FIGURE_DIR

STATION_TERMS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_terms.csv"
OUT_TRANSFER_PREDICTIONS = SUPPLEMENT_DIR / "jshis_station_transfer_stress_predictions.csv"
OUT_TRANSFER_METRICS = SUPPLEMENT_DIR / "jshis_station_transfer_stress_metrics.csv"
OUT_INFLUENTIAL_EVENTS = SUPPLEMENT_DIR / "jshis_influential_event_inventory.csv"
OUT_EVENT_STABILITY = SUPPLEMENT_DIR / "jshis_influential_event_stability.csv"
OUT_AUDIT = SUPPLEMENT_DIR / "jshis_robustness_stress_tests.md"
OUT_FIGURE_PDF = FIGURE_DIR / "supplementary_figure_robustness_stress_tests.pdf"
OUT_FIGURE_PNG = FIGURE_DIR / "supplementary_figure_robustness_stress_tests.png"

REGION_ORDER = [
    "Hokkaido",
    "Northern Honshu",
    "Central-Eastern Japan",
    "Western Honshu-Shikoku",
    "Kyushu-Ryukyu",
]

COMMON_PHYSICAL_FEATURES = [
    "log_vs10",
    "log_avs30",
    "log_d1100",
    "log_d1400",
    "log_d1700",
    "log_d2100",
    "log_dbase",
    "elevation",
    "sensor_depth_glminus",
]


def build_common_feature_model(seed: int) -> Pipeline:
    features = COMMON_PHYSICAL_FEATURES + core.SPATIAL_FEATURES
    preprocess = ColumnTransformer(
        [
            (
                "numeric",
                SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
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


def weighted_correlation(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    weights = np.asarray(weights, dtype=float)
    x_mean = np.average(x, weights=weights)
    y_mean = np.average(y, weights=weights)
    covariance = np.average((x - x_mean) * (y - y_mean), weights=weights)
    x_var = np.average((x - x_mean) ** 2, weights=weights)
    y_var = np.average((y - y_mean) ** 2, weights=weights)
    if x_var <= 0 or y_var <= 0:
        return np.nan
    return float(covariance / np.sqrt(x_var * y_var))


def weighted_quantile(values: np.ndarray, quantile: float, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    positions = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
    return float(np.interp(quantile, positions, values, left=values[0], right=values[-1]))


def assign_macroregion(frame: pd.DataFrame) -> pd.Series:
    lon = pd.to_numeric(frame["lon"], errors="coerce")
    lat = pd.to_numeric(frame["lat"], errors="coerce")
    region = pd.Series("Western Honshu-Shikoku", index=frame.index, dtype="object")
    region.loc[lat.ge(41.3)] = "Hokkaido"
    region.loc[lat.lt(41.3) & lat.ge(37.5)] = "Northern Honshu"
    region.loc[lat.lt(37.5) & lon.ge(137.0)] = "Central-Eastern Japan"
    region.loc[lat.lt(37.5) & lon.lt(132.5)] = "Kyushu-Ryukyu"
    return region


def score_predictions(
    y: np.ndarray,
    prediction: np.ndarray,
    weights: np.ndarray,
) -> dict[str, float]:
    rmse = core.weighted_rmse(y, prediction, weights)
    mae = core.weighted_mae(y, prediction, weights)
    zero_rmse = core.weighted_rmse(y, np.zeros_like(y), weights)
    zero_mae = core.weighted_mae(y, np.zeros_like(y), weights)
    return {
        "weighted_rmse_log10": rmse,
        "weighted_mae_log10": mae,
        "zero_baseline_rmse_log10": zero_rmse,
        "zero_baseline_mae_log10": zero_mae,
        "rmse_reduction_vs_zero_pct": 100.0 * (zero_rmse - rmse) / zero_rmse,
        "mae_reduction_vs_zero_pct": 100.0 * (zero_mae - mae) / zero_mae,
        "weighted_observed_predicted_correlation": weighted_correlation(y, prediction, weights),
        "observed_predicted_spearman": float(pd.Series(y).corr(pd.Series(prediction), method="spearman")),
    }


def fit_and_score_transfer(
    station_terms: pd.DataFrame,
    min_records: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = station_terms[
        station_terms["model"].eq("mf2013_site")
        & station_terms["n_records"].ge(min_records)
        & station_terms["lon"].notna()
        & station_terms["lat"].notna()
    ].copy()
    primary["macroregion"] = assign_macroregion(primary)
    prediction_blocks: list[pd.DataFrame] = []
    metric_rows: list[dict[str, float | int | str]] = []

    for spec in core.PERIODS:
        frame = primary[primary["period_s"].eq(spec.period_s)].copy().reset_index(drop=True)
        y = frame["station_effect_log10"].to_numpy(float)
        weights = frame["n_records"].to_numpy(float)
        validation_groups = {
            "network_transfer": ["K-NET", "KiK-net"],
            "macroregion_leaveout": REGION_ORDER,
        }
        for validation, groups in validation_groups.items():
            group_values = frame["network_label"] if validation == "network_transfer" else frame["macroregion"]
            for held_group in groups:
                test = group_values.eq(held_group).to_numpy()
                train = ~test
                if test.sum() < 30 or train.sum() < 100:
                    continue
                model_specs = [
                    ("physical_hgb", core.build_model(include_space=False, seed=seed)),
                    ("physical_spatial_hgb", core.build_model(include_space=True, seed=seed)),
                    ("common_feature_spatial_hgb", build_common_feature_model(seed=seed)),
                ]
                for model_name, model in model_specs:
                    model.fit(frame.loc[train], y[train], model__sample_weight=weights[train])
                    prediction = model.predict(frame.loc[test])
                    score = score_predictions(y[test], prediction, weights[test])
                    train_group = "other_network" if validation == "network_transfer" else "other_macroregions"
                    metric_rows.append(
                        {
                            "period_s": spec.period_s,
                            "period_code": spec.code,
                            "validation": validation,
                            "model": model_name,
                            "train_group": train_group,
                            "held_group": held_group,
                            "n_train_stations": int(train.sum()),
                            "n_test_stations": int(test.sum()),
                            "train_record_weight": int(weights[train].sum()),
                            "test_record_weight": int(weights[test].sum()),
                            **score,
                        }
                    )
                    block = frame.loc[
                        test,
                        [
                            "siteid2",
                            "site_code",
                            "network_label",
                            "macroregion",
                            "lon",
                            "lat",
                            "n_records",
                            "station_effect_log10",
                        ],
                    ].copy()
                    block["period_s"] = spec.period_s
                    block["period_code"] = spec.code
                    block["validation"] = validation
                    block["model"] = model_name
                    block["train_group"] = train_group
                    block["held_group"] = held_group
                    block["prediction_log10"] = prediction
                    block["prediction_error_log10"] = block["station_effect_log10"] - prediction
                    prediction_blocks.append(block)
        print(f"transfer tests: period={spec.period_s:g}s", flush=True)

    predictions = pd.concat(prediction_blocks, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    return predictions, metrics


def compare_station_terms(
    full: pd.DataFrame,
    perturbed: pd.DataFrame,
    period_s: float,
    perturbation: str,
    removed_event_id: int | None,
    removed_event_count: int,
    min_records: int,
) -> dict[str, float | int | str]:
    full_supported = full[full["n_records"].ge(min_records)][
        ["siteid2", "n_records", "station_effect_log10"]
    ].rename(
        columns={
            "n_records": "full_n_records",
            "station_effect_log10": "full_station_effect_log10",
        }
    )
    perturbed = perturbed.rename(
        columns={
            "n_records": "perturbed_n_records",
            "station_effect_log10": "perturbed_station_effect_log10",
        }
    )
    paired = full_supported.merge(perturbed, on="siteid2", how="inner", validate="one_to_one")
    weights = paired["full_n_records"].to_numpy(float)
    full_values = paired["full_station_effect_log10"].to_numpy(float)
    perturbed_values = paired["perturbed_station_effect_log10"].to_numpy(float)
    reference_shift = float(np.average(perturbed_values - full_values, weights=weights))
    aligned = perturbed_values - reference_shift
    delta = aligned - full_values
    return {
        "period_s": period_s,
        "perturbation": perturbation,
        "removed_event_id": removed_event_id if removed_event_id is not None else -1,
        "removed_event_count": removed_event_count,
        "n_paired_stations": len(paired),
        "record_weight": int(weights.sum()),
        "reference_alignment_shift_log10": reference_shift,
        "weighted_station_term_correlation": weighted_correlation(full_values, aligned, weights),
        "station_term_spearman": float(pd.Series(full_values).corr(pd.Series(aligned), method="spearman")),
        "weighted_mean_absolute_change_log10": float(np.average(np.abs(delta), weights=weights)),
        "weighted_q95_absolute_change_log10": weighted_quantile(np.abs(delta), 0.95, weights),
        "maximum_absolute_change_log10": float(np.max(np.abs(delta))),
    }


def influential_event_stability(
    records: dict[float, pd.DataFrame],
    station_terms: pd.DataFrame,
    source: pd.DataFrame,
    top_n: int,
    min_records: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sa3_records = records[3.0]
    event_counts = (
        sa3_records.groupby("eq_source_id", as_index=False)
        .size()
        .rename(columns={"size": "sa3_record_count"})
        .sort_values(["sa3_record_count", "eq_source_id"], ascending=[False, True])
        .head(top_n)
    )
    source_keep = source[
        ["eq_source_id", "mw", "jem_lat", "jem_lon", "jem_depth", "eq_location_type_id"]
    ].drop_duplicates("eq_source_id")
    inventory = event_counts.merge(source_keep, on="eq_source_id", how="left", validate="one_to_one")
    inventory.insert(0, "influence_rank", np.arange(1, len(inventory) + 1))
    top_ids = set(inventory["eq_source_id"].astype(int))

    rows: list[dict[str, float | int | str]] = []
    full_sa3 = station_terms[
        station_terms["model"].eq("mf2013_site") & np.isclose(station_terms["period_s"], 3.0)
    ]
    for event_id in inventory["eq_source_id"].astype(int):
        perturbed, _, _ = core.two_way_decomposition(
            sa3_records[~sa3_records["eq_source_id"].eq(event_id)],
            "residual_site",
        )
        rows.append(
            compare_station_terms(
                full_sa3,
                perturbed,
                period_s=3.0,
                perturbation="single_event_removal",
                removed_event_id=event_id,
                removed_event_count=1,
                min_records=min_records,
            )
        )
        print(f"event influence: removed={event_id}", flush=True)

    for spec in core.PERIODS:
        full = station_terms[
            station_terms["model"].eq("mf2013_site") & station_terms["period_s"].eq(spec.period_s)
        ]
        period_records = records[spec.period_s]
        perturbed, _, _ = core.two_way_decomposition(
            period_records[~period_records["eq_source_id"].isin(top_ids)],
            "residual_site",
        )
        rows.append(
            compare_station_terms(
                full,
                perturbed,
                period_s=spec.period_s,
                perturbation=f"joint_top{top_n}_removal",
                removed_event_id=None,
                removed_event_count=top_n,
                min_records=min_records,
            )
        )
        print(f"event influence: period={spec.period_s:g}s joint top {top_n}", flush=True)
    return inventory, pd.DataFrame(rows)


def plot_figure(
    transfer_metrics: pd.DataFrame,
    event_inventory: pd.DataFrame,
    event_stability: pd.DataFrame,
    top_n: int,
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
    colors = {"K-NET": "#126782", "KiK-net": "#F28E2B", "joint": "#7B2CBF"}
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2), constrained_layout=True)

    network = transfer_metrics[
        transfer_metrics["validation"].eq("network_transfer")
        & transfer_metrics["model"].eq("common_feature_spatial_hgb")
    ]
    ax = axes[0, 0]
    for held_group, block in network.groupby("held_group"):
        block = block.sort_values("period_s")
        ax.semilogx(
            block["period_s"],
            block["weighted_observed_predicted_correlation"],
            marker="o",
            color=colors[held_group],
            label=f"Held {held_group}",
        )
    ax.set(xlabel="Oscillator period (s)", ylabel="Weighted correlation", title="a  Cross-network transfer")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.grid(alpha=0.2, linewidth=0.5)
    ax.legend(frameon=False)

    region = transfer_metrics[
        transfer_metrics["validation"].eq("macroregion_leaveout")
        & transfer_metrics["model"].eq("common_feature_spatial_hgb")
    ]
    matrix = region.pivot(index="held_group", columns="period_s", values="rmse_reduction_vs_zero_pct")
    matrix = matrix.reindex(index=REGION_ORDER, columns=[spec.period_s for spec in core.PERIODS])
    ax = axes[0, 1]
    vmax = max(40.0, float(np.nanmax(np.abs(matrix.to_numpy()))))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix.iloc[row, col]
            ax.text(col, row, f"{value:.0f}", ha="center", va="center", fontsize=6.7)
    ax.set_xticks(range(matrix.shape[1]), [f"{value:g}" for value in matrix.columns])
    ax.set_yticks(range(matrix.shape[0]), matrix.index)
    ax.set(xlabel="Oscillator period (s)", title="b  Leave-macroregion-out RMSE gain (%)")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.047, pad=0.02)
    colorbar.set_label("Gain vs zero baseline (%)")

    single = event_stability[event_stability["perturbation"].eq("single_event_removal")].copy()
    single = single.merge(
        event_inventory[["influence_rank", "eq_source_id", "sa3_record_count"]],
        left_on="removed_event_id",
        right_on="eq_source_id",
        how="left",
        validate="one_to_one",
    ).sort_values("influence_rank")
    ax = axes[1, 0]
    ax.plot(
        single["influence_rank"],
        single["weighted_station_term_correlation"],
        marker="o",
        color="#126782",
    )
    joint_sa3 = event_stability[
        event_stability["perturbation"].eq(f"joint_top{top_n}_removal")
        & np.isclose(event_stability["period_s"], 3.0)
    ].iloc[0]
    ax.axhline(
        joint_sa3["weighted_station_term_correlation"],
        color=colors["joint"],
        linestyle="--",
        label=f"Joint top-{top_n} removal",
    )
    ax.set(
        xlabel="Event rank by SA(3.0 s) record count",
        ylabel="Station-term correlation",
        title="c  SA(3.0 s) influential-event deletion",
    )
    lower = min(single["weighted_station_term_correlation"].min(), joint_sa3["weighted_station_term_correlation"])
    ax.set_ylim(max(0.9, lower - 0.005), 1.0005)
    ax.grid(alpha=0.2, linewidth=0.5)
    ax.legend(frameon=False, loc="lower right")

    joint = event_stability[event_stability["perturbation"].eq(f"joint_top{top_n}_removal")].sort_values("period_s")
    ax = axes[1, 1]
    line1 = ax.semilogx(
        joint["period_s"],
        joint["weighted_station_term_correlation"],
        marker="o",
        color=colors["joint"],
        label="Weighted correlation",
    )[0]
    ax.set(xlabel="Oscillator period (s)", ylabel="Weighted correlation", title=f"d  Joint top-{top_n} event removal")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.grid(alpha=0.2, linewidth=0.5)
    twin = ax.twinx()
    line2 = twin.semilogx(
        joint["period_s"],
        joint["weighted_q95_absolute_change_log10"],
        marker="s",
        color="#F28E2B",
        label="95th-percentile change",
    )[0]
    twin.set_ylabel("95th-percentile absolute change (log10)")
    ax.legend([line1, line2], [line1.get_label(), line2.get_label()], frameon=False, loc="best")

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    transfer_metrics: pd.DataFrame,
    event_stability: pd.DataFrame,
    top_n: int,
) -> None:
    network = transfer_metrics[
        transfer_metrics["validation"].eq("network_transfer")
        & transfer_metrics["model"].eq("physical_spatial_hgb")
    ]
    common_network = transfer_metrics[
        transfer_metrics["validation"].eq("network_transfer")
        & transfer_metrics["model"].eq("common_feature_spatial_hgb")
    ]
    regions = transfer_metrics[
        transfer_metrics["validation"].eq("macroregion_leaveout")
        & transfer_metrics["model"].eq("physical_spatial_hgb")
    ]
    common_regions = transfer_metrics[
        transfer_metrics["validation"].eq("macroregion_leaveout")
        & transfer_metrics["model"].eq("common_feature_spatial_hgb")
    ]
    sa3_network = network[np.isclose(network["period_s"], 3.0)]
    sa3_common_network = common_network[np.isclose(common_network["period_s"], 3.0)]
    sa3_region = regions[np.isclose(regions["period_s"], 3.0)]
    sa3_common_region = common_regions[np.isclose(common_regions["period_s"], 3.0)]
    single = event_stability[event_stability["perturbation"].eq("single_event_removal")]
    joint = event_stability[event_stability["perturbation"].eq(f"joint_top{top_n}_removal")]
    sa3_joint = joint[np.isclose(joint["period_s"], 3.0)].iloc[0]
    lines = [
        "# Network, geography, and influential-event stress tests",
        "",
        "## Validation design",
        "",
        "- Cross-network transfer trains the station model on one network and evaluates it on the other without using target-network station terms.",
        "- Fixed macroregion tests leave out one of five coordinate-defined regions. These labels are reproducible geographic partitions, not administrative regions.",
        "- The common-feature sensitivity excludes VS20, VS30 and network identity because their completeness differs between K-NET and KiK-net; it retains VS10, AVS30, sedimentary depths, elevation, sensor depth and geographic variables.",
        f"- Influential-event tests remove each of the {top_n} events with the largest SA(3.0 s) record counts and then remove all {top_n} jointly.",
        "- Event-deletion comparisons align the perturbed and full station-term references by their record-weighted mean difference before computing changes.",
        "",
        "## Results",
        "",
        f"- SA(3.0 s) cross-network weighted correlations span {sa3_network['weighted_observed_predicted_correlation'].min():.3f}--{sa3_network['weighted_observed_predicted_correlation'].max():.3f}.",
        f"- SA(3.0 s) cross-network RMSE gains span {sa3_network['rmse_reduction_vs_zero_pct'].min():.1f}%--{sa3_network['rmse_reduction_vs_zero_pct'].max():.1f}%.",
        f"- With common-coverage features, SA(3.0 s) cross-network correlations span {sa3_common_network['weighted_observed_predicted_correlation'].min():.3f}--{sa3_common_network['weighted_observed_predicted_correlation'].max():.3f}, and RMSE gains span {sa3_common_network['rmse_reduction_vs_zero_pct'].min():.1f}%--{sa3_common_network['rmse_reduction_vs_zero_pct'].max():.1f}%.",
        f"- SA(3.0 s) leave-macroregion-out RMSE gains span {sa3_region['rmse_reduction_vs_zero_pct'].min():.1f}%--{sa3_region['rmse_reduction_vs_zero_pct'].max():.1f}%.",
        f"- Common-feature leave-macroregion-out SA(3.0 s) RMSE gains span {sa3_common_region['rmse_reduction_vs_zero_pct'].min():.1f}%--{sa3_common_region['rmse_reduction_vs_zero_pct'].max():.1f}%.",
        f"- Individual high-count event deletion retains SA(3.0 s) station-term correlations of {single['weighted_station_term_correlation'].min():.4f}--{single['weighted_station_term_correlation'].max():.4f}.",
        f"- Joint top-{top_n} deletion retains an SA(3.0 s) correlation of {sa3_joint['weighted_station_term_correlation']:.4f}, with a weighted 95th-percentile absolute change of {sa3_joint['weighted_q95_absolute_change_log10']:.3f} log10 units.",
        f"- Joint top-{top_n} deletion correlations across all periods span {joint['weighted_station_term_correlation'].min():.4f}--{joint['weighted_station_term_correlation'].max():.4f}.",
        "",
        "## Interpretation boundary",
        "",
        "The transfer tests reuse station terms estimated from the same observational archive; they are model-generalization tests rather than independent-network replications. Event deletion tests sensitivity to high-count events and do not replace a fully independent earthquake catalogue.",
        "",
    ]
    OUT_AUDIT.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=core.DEFAULT_COEFFICIENTS)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--min-station-records", type=int, default=20)
    parser.add_argument("--top-events", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    station_terms = pd.read_csv(STATION_TERMS, low_memory=False)
    transfer_predictions, transfer_metrics = fit_and_score_transfer(
        station_terms,
        min_records=args.min_station_records,
        seed=args.seed,
    )
    coefficients = core.load_coefficients(args.coefficients)
    site, source = core.load_metadata(args.flatfile)
    records = core.load_record_residuals(args.flatfile, coefficients, site, source, args.chunksize)
    event_inventory, event_stability = influential_event_stability(
        records,
        station_terms,
        source,
        top_n=args.top_events,
        min_records=args.min_station_records,
    )

    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    transfer_predictions.to_csv(OUT_TRANSFER_PREDICTIONS, index=False)
    transfer_metrics.to_csv(OUT_TRANSFER_METRICS, index=False)
    event_inventory.to_csv(OUT_INFLUENTIAL_EVENTS, index=False)
    event_stability.to_csv(OUT_EVENT_STABILITY, index=False)
    plot_figure(transfer_metrics, event_inventory, event_stability, args.top_events)
    write_audit(transfer_metrics, event_stability, args.top_events)
    print(f"wrote {OUT_AUDIT}", flush=True)


if __name__ == "__main__":
    main()
