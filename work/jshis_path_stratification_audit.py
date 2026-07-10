#!/usr/bin/env python3
"""Audit MF2013 station terms across path and source strata."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsmr

import jshis_event_adjusted_station_model as core
from jshis_robustness_stress_tests import weighted_correlation, weighted_quantile


OUT_SUMMARY = core.SUPPLEMENT_DIR / "jshis_path_stratification_stability.csv"
OUT_STATIONS = core.SUPPLEMENT_DIR / "jshis_path_stratification_station_terms.csv"
OUT_REPEATABILITY = core.SUPPLEMENT_DIR / "jshis_path_stratification_split_half_repeatability.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_path_stratification_audit.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "figure_path_stratification.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "figure_path_stratification.png"

AZIMUTH_GROUPS = ["000-090 deg", "090-180 deg", "180-270 deg", "270-360 deg"]
DISTANCE_GROUPS = ["0-50 km", "50-100 km", "100-200 km", "200-300 km"]
SOURCE_GROUPS = ["Crustal", "Interplate", "Intraplate"]
GROUPS = {
    "hypocentral_azimuth": AZIMUTH_GROUPS,
    "fault_distance": DISTANCE_GROUPS,
    "source_type": SOURCE_GROUPS,
}


def hypocentral_azimuth(frame: pd.DataFrame) -> np.ndarray:
    source_lat = np.radians(frame["source_lat"].to_numpy(float))
    site_lat = np.radians(frame["site_lat"].to_numpy(float))
    delta_lon = np.radians(frame["site_lon"].to_numpy(float) - frame["source_lon"].to_numpy(float))
    x = np.sin(delta_lon) * np.cos(site_lat)
    y = np.cos(source_lat) * np.sin(site_lat) - np.sin(source_lat) * np.cos(site_lat) * np.cos(delta_lon)
    return (np.degrees(np.arctan2(x, y)) + 360.0) % 360.0


def add_record_strata(
    records: pd.DataFrame,
    site_lookup: pd.DataFrame,
    source_lookup: pd.DataFrame,
) -> pd.DataFrame:
    frame = records.merge(site_lookup, on="siteid2", how="left", validate="many_to_one")
    frame = frame.merge(source_lookup, on="eq_source_id", how="left", validate="many_to_one")
    frame["hypocentral_azimuth_deg"] = hypocentral_azimuth(frame)
    frame["hypocentral_azimuth"] = pd.cut(
        frame["hypocentral_azimuth_deg"],
        bins=[0.0, 90.0, 180.0, 270.0, 360.0],
        labels=AZIMUTH_GROUPS,
        right=False,
        include_lowest=True,
    ).astype("object")
    frame["fault_distance"] = pd.cut(
        pd.to_numeric(frame["fault_dist"], errors="coerce"),
        bins=[0.0, 50.0, 100.0, 200.0, np.inf],
        labels=DISTANCE_GROUPS,
        right=False,
        include_lowest=True,
    ).astype("object")
    frame["source_type"] = frame["eq_location_type_id"].map(
        {1: "Crustal", 2: "Interplate", 3: "Intraplate"}
    )
    return frame


def solve_connected_component(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float | int]]:
    clean = frame[["siteid2", "eq_source_id", "residual_site"]].dropna().copy()
    station_codes, station_ids = pd.factorize(clean["siteid2"], sort=True)
    event_codes, event_ids = pd.factorize(clean["eq_source_id"], sort=True)
    n_records = len(clean)
    n_stations = len(station_ids)
    n_events = len(event_ids)
    rows = [np.arange(n_records)]
    columns = [np.zeros(n_records, dtype=int)]
    data = [np.ones(n_records)]
    station_mask = station_codes > 0
    event_mask = event_codes > 0
    rows.extend([np.flatnonzero(station_mask), np.flatnonzero(event_mask)])
    columns.extend([station_codes[station_mask], n_stations + event_codes[event_mask] - 1])
    data.extend([np.ones(station_mask.sum()), np.ones(event_mask.sum())])
    design = coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(columns))),
        shape=(n_records, n_stations + n_events - 1),
    ).tocsr()
    values = clean["residual_site"].to_numpy(float)
    result = lsmr(
        design,
        values,
        atol=1.0e-12,
        btol=1.0e-12,
        maxiter=20_000,
    )
    solution, stop_code, iterations, _, normal_residual, _, condition = result[:7]
    station_effect = np.zeros(n_stations)
    event_effect = np.zeros(n_events)
    station_effect[1:] = solution[1:n_stations]
    event_effect[1:] = solution[n_stations:]
    station_counts = np.bincount(station_codes).astype(float)
    event_counts = np.bincount(event_codes).astype(float)
    station_effect -= np.average(station_effect, weights=station_counts)
    event_effect -= np.average(event_effect, weights=event_counts)
    terms = pd.DataFrame(
        {
            "siteid2": pd.Series(station_ids, dtype="Int64"),
            "n_records": station_counts.astype(int),
            "station_effect_log10": station_effect,
        }
    )
    relative_normal_residual = normal_residual / max(np.linalg.norm(values), np.finfo(float).eps)
    return terms, {
        "solver_stop_code": int(stop_code),
        "solver_iterations": int(iterations),
        "solver_condition_estimate": float(condition),
        "solver_relative_normal_residual": float(relative_normal_residual),
    }


def component_station_terms(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    clean = frame[["siteid2", "eq_source_id", "residual_site"]].dropna().copy()
    station_codes, station_ids = pd.factorize(clean["siteid2"], sort=True)
    event_codes, event_ids = pd.factorize(clean["eq_source_id"], sort=True)
    n_stations = len(station_ids)
    rows = np.concatenate([station_codes, n_stations + event_codes])
    columns = np.concatenate([n_stations + event_codes, station_codes])
    graph = coo_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, columns)),
        shape=(n_stations + len(event_ids), n_stations + len(event_ids)),
    ).tocsr()
    n_components, labels = connected_components(graph, directed=False)
    record_components = labels[station_codes]
    node_counts = np.bincount(labels)
    term_blocks = []
    metric_rows = []
    for component in range(n_components):
        block = clean.loc[record_components == component]
        terms, metrics = solve_connected_component(block)
        terms["graph_component"] = component
        terms["component_node_count"] = int(node_counts[component])
        term_blocks.append(terms)
        metric_rows.append(
            {
                "graph_component": component,
                "component_node_count": int(node_counts[component]),
                "component_records": len(block),
                **metrics,
            }
        )
    return pd.concat(term_blocks, ignore_index=True), pd.DataFrame(metric_rows)


def align_component_terms(
    stratum_terms: pd.DataFrame,
    full_terms: pd.DataFrame,
    min_stratum_records: int,
    min_component_stations: int,
) -> pd.DataFrame | None:
    supported = stratum_terms[stratum_terms["n_records"].ge(min_stratum_records)].rename(
        columns={
            "n_records": "stratum_n_records",
            "station_effect_log10": "stratum_station_effect_log10",
        }
    )
    paired_blocks = []
    for _, component_terms in supported.groupby("graph_component"):
        paired = full_terms.merge(component_terms, on="siteid2", how="inner", validate="one_to_one")
        if len(paired) < min_component_stations:
            continue
        weights = paired["stratum_n_records"].to_numpy(float)
        full_values = paired["full_station_effect_log10"].to_numpy(float)
        stratum_values = paired["stratum_station_effect_log10"].to_numpy(float)
        reference_shift = float(np.average(stratum_values - full_values, weights=weights))
        paired["reference_alignment_shift_log10"] = reference_shift
        paired["aligned_stratum_station_effect_log10"] = stratum_values - reference_shift
        paired["aligned_difference_log10"] = paired["aligned_stratum_station_effect_log10"] - full_values
        paired_blocks.append(paired)
    return pd.concat(paired_blocks, ignore_index=True) if paired_blocks else None


def weighted_group_centered_correlation(
    first: np.ndarray,
    second: np.ndarray,
    weights: np.ndarray,
    groups: pd.Series,
) -> float:
    """Correlate values after removing an independent offset in each graph block."""
    first_centered = np.empty_like(first, dtype=float)
    second_centered = np.empty_like(second, dtype=float)
    group_values = groups.to_numpy()
    for group in pd.unique(group_values):
        selected = group_values == group
        group_weights = weights[selected]
        first_centered[selected] = first[selected] - np.average(first[selected], weights=group_weights)
        second_centered[selected] = second[selected] - np.average(second[selected], weights=group_weights)
    return weighted_correlation(first_centered, second_centered, weights)


def compare_stratum(
    subset: pd.DataFrame,
    full_terms: pd.DataFrame,
    period_s: float,
    period_code: str,
    dimension: str,
    group: str,
    min_stratum_records: int,
    min_component_stations: int,
) -> tuple[dict[str, float | int | str], pd.DataFrame] | None:
    stratum_terms, component_metrics = component_station_terms(subset)
    paired = align_component_terms(
        stratum_terms,
        full_terms,
        min_stratum_records=min_stratum_records,
        min_component_stations=min_component_stations,
    )
    if paired is None:
        return None
    weights = paired["stratum_n_records"].to_numpy(float)
    full_values = paired["full_station_effect_log10"].to_numpy(float)
    aligned = paired["aligned_stratum_station_effect_log10"].to_numpy(float)
    differences = aligned - full_values

    paired["period_s"] = period_s
    paired["period_code"] = period_code
    paired["stratification"] = dimension
    paired["stratum"] = group

    metrics: dict[str, float | int | str] = {
        "period_s": period_s,
        "period_code": period_code,
        "stratification": dimension,
        "stratum": group,
        "n_records": len(subset),
        "n_events": subset["eq_source_id"].nunique(),
        "n_stations_raw": subset["siteid2"].nunique(),
        "n_paired_stations": len(paired),
        "n_graph_components": len(component_metrics),
        "n_components_retained": paired["graph_component"].nunique(),
        "largest_component_node_fraction": (
            component_metrics["component_node_count"].max() / component_metrics["component_node_count"].sum()
        ),
        "paired_record_weight": int(weights.sum()),
        "median_stratum_records_per_station": float(np.median(weights)),
        "weighted_reference_alignment_shift_log10": float(
            np.average(paired["reference_alignment_shift_log10"], weights=weights)
        ),
        "weighted_station_term_correlation": weighted_correlation(full_values, aligned, weights),
        "weighted_within_component_correlation": weighted_group_centered_correlation(
            full_values,
            aligned,
            weights,
            paired["graph_component"],
        ),
        "station_term_spearman": float(pd.Series(full_values).corr(pd.Series(aligned), method="spearman")),
        "weighted_rmse_difference_log10": core.weighted_rmse(full_values, aligned, weights),
        "weighted_mean_absolute_difference_log10": core.weighted_mae(full_values, aligned, weights),
        "weighted_q95_absolute_difference_log10": weighted_quantile(np.abs(differences), 0.95, weights),
        "maximum_absolute_difference_log10": float(np.max(np.abs(differences))),
        "solver_max_iterations": int(component_metrics["solver_iterations"].max()),
        "solver_max_condition_estimate": float(component_metrics["solver_condition_estimate"].max()),
        "solver_nonconverged_components": int(
            (~component_metrics["solver_stop_code"].isin([1, 2, 4, 5])).sum()
        ),
        "solver_max_relative_normal_residual": float(
            component_metrics["solver_relative_normal_residual"].max()
        ),
    }
    return metrics, paired


def evaluate(
    records: dict[float, pd.DataFrame],
    station_terms: pd.DataFrame,
    site: pd.DataFrame,
    source: pd.DataFrame,
    min_full_records: int,
    min_stratum_records: int,
    min_component_stations: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    site_lookup = site[["siteid2", "lon", "lat"]].drop_duplicates("siteid2").rename(
        columns={"lon": "site_lon", "lat": "site_lat"}
    )
    source_lookup = source[
        ["eq_source_id", "jem_lon", "jem_lat", "eq_location_type_id"]
    ].drop_duplicates("eq_source_id").rename(
        columns={"jem_lon": "source_lon", "jem_lat": "source_lat"}
    )
    summary_rows: list[dict[str, float | int | str]] = []
    station_rows: list[pd.DataFrame] = []

    for spec in core.PERIODS:
        frame = add_record_strata(records[spec.period_s], site_lookup, source_lookup)
        full = station_terms[
            station_terms["period_s"].eq(spec.period_s)
            & station_terms["model"].eq("mf2013_site")
            & station_terms["n_records"].ge(min_full_records)
        ][["siteid2", "n_records", "station_effect_log10"]].rename(
            columns={
                "n_records": "full_n_records",
                "station_effect_log10": "full_station_effect_log10",
            }
        )
        for dimension, groups in GROUPS.items():
            for group in groups:
                subset = frame[frame[dimension].eq(group)]
                result = compare_stratum(
                    subset,
                    full,
                    period_s=spec.period_s,
                    period_code=spec.code,
                    dimension=dimension,
                    group=group,
                    min_stratum_records=min_stratum_records,
                    min_component_stations=min_component_stations,
                )
                if result is None:
                    continue
                metrics, paired = result
                summary_rows.append(metrics)
                station_rows.append(paired)
        print(f"path stratification: period={spec.period_s:g}s", flush=True)

    return pd.DataFrame(summary_rows), pd.concat(station_rows, ignore_index=True)


def split_half_repeatability(
    records: dict[float, pd.DataFrame],
    station_terms: pd.DataFrame,
    site: pd.DataFrame,
    source: pd.DataFrame,
    min_full_records: int,
    min_half_records: int,
    min_component_stations: int,
    n_splits: int,
    seed: int,
) -> pd.DataFrame:
    site_lookup = site[["siteid2", "lon", "lat"]].drop_duplicates("siteid2").rename(
        columns={"lon": "site_lon", "lat": "site_lat"}
    )
    source_lookup = source[
        ["eq_source_id", "jem_lon", "jem_lat", "eq_location_type_id"]
    ].drop_duplicates("eq_source_id").rename(
        columns={"jem_lon": "source_lon", "jem_lat": "source_lat"}
    )
    rows: list[dict[str, float | int | str]] = []
    for spec in [item for item in core.PERIODS if item.period_s in (1.0, 2.0, 3.0)]:
        frame = add_record_strata(records[spec.period_s], site_lookup, source_lookup)
        full = station_terms[
            station_terms["period_s"].eq(spec.period_s)
            & station_terms["model"].eq("mf2013_site")
            & station_terms["n_records"].ge(min_full_records)
        ][["siteid2", "n_records", "station_effect_log10"]].rename(
            columns={
                "n_records": "full_n_records",
                "station_effect_log10": "full_station_effect_log10",
            }
        )
        for dimension_index, (dimension, groups) in enumerate(GROUPS.items()):
            for group_index, group in enumerate(groups):
                subset = frame[frame[dimension].eq(group)]
                events = np.sort(subset["eq_source_id"].dropna().unique())
                for split in range(n_splits):
                    rng = np.random.default_rng(
                        seed
                        + int(round(spec.period_s * 1000)) * 1000
                        + dimension_index * 100
                        + group_index * 10
                        + split
                    )
                    shuffled = rng.permutation(events)
                    event_halves = (shuffled[::2], shuffled[1::2])
                    aligned_halves = []
                    for event_half in event_halves:
                        terms, _ = component_station_terms(subset[subset["eq_source_id"].isin(event_half)])
                        aligned = align_component_terms(
                            terms,
                            full,
                            min_stratum_records=min_half_records,
                            min_component_stations=min_component_stations,
                        )
                        if aligned is None:
                            break
                        aligned_halves.append(aligned)
                    if len(aligned_halves) != 2:
                        continue
                    half_a = aligned_halves[0][
                        [
                            "siteid2",
                            "stratum_n_records",
                            "full_station_effect_log10",
                            "aligned_stratum_station_effect_log10",
                            "graph_component",
                        ]
                    ].rename(
                        columns={
                            "stratum_n_records": "half_a_n_records",
                            "aligned_stratum_station_effect_log10": "half_a_station_effect_log10",
                            "graph_component": "half_a_graph_component",
                        }
                    )
                    half_b = aligned_halves[1][
                        [
                            "siteid2",
                            "stratum_n_records",
                            "aligned_stratum_station_effect_log10",
                            "graph_component",
                        ]
                    ].rename(
                        columns={
                            "stratum_n_records": "half_b_n_records",
                            "aligned_stratum_station_effect_log10": "half_b_station_effect_log10",
                            "graph_component": "half_b_graph_component",
                        }
                    )
                    paired = half_a.merge(half_b, on="siteid2", how="inner", validate="one_to_one")
                    if len(paired) < 30:
                        continue
                    weights = paired[["half_a_n_records", "half_b_n_records"]].min(axis=1).to_numpy(float)
                    full_values = paired["full_station_effect_log10"].to_numpy(float)
                    half_a_values = paired["half_a_station_effect_log10"].to_numpy(float)
                    half_b_values = paired["half_b_station_effect_log10"].to_numpy(float)
                    component_pairs = (
                        paired["half_a_graph_component"].astype(str)
                        + ":"
                        + paired["half_b_graph_component"].astype(str)
                    )
                    rows.append(
                        {
                            "period_s": spec.period_s,
                            "period_code": spec.code,
                            "stratification": dimension,
                            "stratum": group,
                            "split": split,
                            "n_events_half_a": len(event_halves[0]),
                            "n_events_half_b": len(event_halves[1]),
                            "n_paired_stations": len(paired),
                            "n_components_half_a_retained": paired["half_a_graph_component"].nunique(),
                            "n_components_half_b_retained": paired["half_b_graph_component"].nunique(),
                            "n_component_pairs": component_pairs.nunique(),
                            "paired_record_weight": int(weights.sum()),
                            "weighted_split_half_correlation": weighted_correlation(
                                half_a_values, half_b_values, weights
                            ),
                            "weighted_split_half_correlation_within_component": weighted_group_centered_correlation(
                                half_a_values,
                                half_b_values,
                                weights,
                                component_pairs,
                            ),
                            "split_half_spearman": float(
                                pd.Series(half_a_values).corr(pd.Series(half_b_values), method="spearman")
                            ),
                            "weighted_split_half_rmse_log10": core.weighted_rmse(
                                half_a_values, half_b_values, weights
                            ),
                            "weighted_half_a_full_correlation": weighted_correlation(
                                half_a_values, full_values, weights
                            ),
                            "weighted_half_b_full_correlation": weighted_correlation(
                                half_b_values, full_values, weights
                            ),
                        }
                    )
        print(f"path split-half: period={spec.period_s:g}s", flush=True)
    return pd.DataFrame(rows)


def save_figure(summary: pd.DataFrame, repeatability: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 6.8,
            "axes.linewidth": 0.7,
            "savefig.dpi": 300,
        }
    )
    titles = {
        "hypocentral_azimuth": "a  Hypocentral-bearing sectors",
        "fault_distance": "b  Shortest-fault-distance bins",
        "source_type": "c  Source classes",
    }
    colors = ["#126782", "#E76F51", "#2A9D8F", "#6A4C93"]
    markers = ["o", "s", "^", "D"]
    repeatability_mean = (
        repeatability.groupby(["period_s", "stratification", "stratum"], as_index=False)
        ["weighted_split_half_correlation_within_component"]
        .mean()
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), sharey=True, constrained_layout=True)
    for ax, (dimension, groups) in zip(axes, GROUPS.items(), strict=True):
        block = summary[summary["stratification"].eq(dimension)]
        for index, group in enumerate(groups):
            line = block[block["stratum"].eq(group)].sort_values("period_s")
            if line.empty:
                continue
            ax.semilogx(
                line["period_s"],
                line["weighted_station_term_correlation"],
                color=colors[index],
                marker=markers[index],
                linewidth=1.1,
                markersize=3.5,
                label=group,
            )
            repeat = repeatability_mean[
                repeatability_mean["stratification"].eq(dimension)
                & repeatability_mean["stratum"].eq(group)
            ].sort_values("period_s")
            ax.semilogx(
                repeat["period_s"],
                repeat["weighted_split_half_correlation_within_component"],
                color=colors[index],
                marker=markers[index],
                linestyle="--",
                linewidth=0.9,
                markersize=3.0,
                alpha=0.75,
            )
        ax.set_title(titles[dimension], loc="left")
        ax.set_xlabel("Oscillator period (s)")
        ax.set_xticks([0.1, 0.2, 0.5, 1.0, 2.0, 5.0], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
        ax.grid(alpha=0.2, linewidth=0.5)
        ax.legend(frameon=False, loc="lower right")
    axes[0].set_ylabel("Weighted station-term correlation")
    minimum = float(summary["weighted_station_term_correlation"].min())
    axes[0].set_ylim(min(-0.05, minimum - 0.05), 1.01)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    summary: pd.DataFrame,
    repeatability: pd.DataFrame,
    min_stratum_records: int,
    min_half_records: int,
    min_component_stations: int,
) -> None:
    focus = summary[summary["period_s"].isin([1.0, 2.0, 3.0])]
    lines = [
        "# Path and source stratification audit",
        "",
        "## Design",
        "",
        "- The primary RotD100 residuals are re-estimated within four hypocentral-bearing sectors, four shortest-fault-distance bins and three J-SHIS source classes.",
        "- Each stratum is separated into connected station-event graph components and solved by sparse least squares. This preserves the primary additive event-station model without imposing offsets between disconnected components.",
        f"- A station requires at least {min_stratum_records} records within a stratum. A connected component requires at least {min_component_stations} supported stations. Each retained component is aligned to the full-sample field by its record-weighted mean difference before comparison.",
        f"- Sampling reliability is evaluated with repeated disjoint event halves at 1, 2 and 3 s. Each station requires at least {min_half_records} records per half. We report both the reference-aligned correlation and a correlation after removing an independent weighted offset from every intersecting component pair.",
        "",
        "## Results at 1-3 s",
        "",
    ]
    for dimension, label in [
        ("hypocentral_azimuth", "Hypocentral-bearing sectors"),
        ("fault_distance", "Shortest-fault-distance bins"),
        ("source_type", "Source classes"),
    ]:
        block = focus[focus["stratification"].eq(dimension)]
        lines.append(
            f"- {label}: weighted station-term correlations span "
            f"{block['weighted_station_term_correlation'].min():.3f}--{block['weighted_station_term_correlation'].max():.3f}; "
            f"the largest weighted 95th-percentile absolute difference is "
            f"{block['weighted_q95_absolute_difference_log10'].max():.3f} log10 units, with at least "
            f"{int(block['n_paired_stations'].min()):,} paired stations."
        )
    repeatability_means = (
        repeatability.groupby(["period_s", "stratification", "stratum"], as_index=False)
        .agg(
            mean_split_half_correlation=(
                "weighted_split_half_correlation_within_component",
                "mean",
            ),
            minimum_paired_stations=("n_paired_stations", "min"),
        )
    )
    lines.extend(["", "## Split-half sampling reliability", ""])
    for dimension, label in [
        ("hypocentral_azimuth", "Hypocentral-bearing sectors"),
        ("fault_distance", "Shortest-fault-distance bins"),
        ("source_type", "Source classes"),
    ]:
        block = repeatability_means[repeatability_means["stratification"].eq(dimension)]
        lines.append(
            f"- {label}: mean within-component split-half correlations span "
            f"{block['mean_split_half_correlation'].min():.3f}--{block['mean_split_half_correlation'].max():.3f}; "
            f"each split retains at least {int(block['minimum_paired_stations'].min()):,} paired stations."
        )
    component_free = (
        repeatability.groupby(["stratification", "stratum"], as_index=False)
        .agg(
            reference_aligned=("weighted_split_half_correlation", "mean"),
            within_component=("weighted_split_half_correlation_within_component", "mean"),
        )
    )
    component_free["alignment_effect"] = (
        component_free["reference_aligned"] - component_free["within_component"]
    ).abs()
    lines.append(
        "- Largest absolute change after removing component-pair offsets: "
        f"{component_free['alignment_effect'].max():.3f}."
    )
    lines.extend(
        [
            f"- Maximum relative normal-equation residual across all component solves: {summary['solver_max_relative_normal_residual'].max():.3e}.",
            f"- Component solves outside accepted LSMR stop codes: {int(summary['solver_nonconverged_components'].sum())}.",
            "",
            "## Interpretation boundary",
            "",
            "The bearing is the initial great-circle direction from the event horizontal coordinates to the station; it is a reproducible directional proxy rather than a finite-fault ray path. Component-specific alignment removes offsets that are not identifiable in disconnected station-event graphs. Reported split-half correlations are invariant to those offsets because each intersecting component pair is centred independently. Distance and source-class subsets contain fewer records per station, so their differences combine sampling uncertainty with any unresolved path dependence. This audit tests whether the station field is dominated by one path sector or source class. It does not identify a path-specific nonergodic term or fully separate site and path effects.",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    core.SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    core.FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    coefficients = core.load_coefficients(args.coefficients)
    site, source = core.load_metadata(args.flatfile)
    records = core.load_record_residuals(args.flatfile, coefficients, site, source, args.chunksize)
    station_terms = pd.read_csv(args.station_terms)
    summary, stations = evaluate(
        records,
        station_terms,
        site,
        source,
        min_full_records=args.min_full_records,
        min_stratum_records=args.min_stratum_records,
        min_component_stations=args.min_component_stations,
    )
    repeatability = split_half_repeatability(
        records,
        station_terms,
        site,
        source,
        min_full_records=args.min_full_records,
        min_half_records=args.min_half_records,
        min_component_stations=args.min_component_stations,
        n_splits=args.repeatability_splits,
        seed=args.seed,
    )
    summary.to_csv(OUT_SUMMARY, index=False)
    stations.to_csv(OUT_STATIONS, index=False)
    repeatability.to_csv(OUT_REPEATABILITY, index=False)
    save_figure(summary, repeatability)
    write_audit(
        summary,
        repeatability,
        args.min_stratum_records,
        args.min_half_records,
        args.min_component_stations,
    )
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=core.DEFAULT_COEFFICIENTS)
    parser.add_argument("--station-terms", type=Path, default=core.OUT_STATION_TERMS)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--min-full-records", type=int, default=20)
    parser.add_argument("--min-stratum-records", type=int, default=10)
    parser.add_argument("--min-half-records", type=int, default=5)
    parser.add_argument("--min-component-stations", type=int, default=10)
    parser.add_argument("--repeatability-splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
