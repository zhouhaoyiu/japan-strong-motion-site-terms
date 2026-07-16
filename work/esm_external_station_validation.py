#!/usr/bin/env python3
"""Independent 1--3 s station-field replication with ESM v3 records."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openquake.hazardlib.gsim.bindi_2014 import BindiEtAl2014RhypEC8NoSOF
from openquake.hazardlib.gsim.kotha_2020 import KothaEtAl2020Site
from openquake.hazardlib.imt import SA
from sklearn.cluster import KMeans
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

import jshis_event_adjusted_station_model as core
from jshis_path_stratification_audit import component_station_terms


ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = ROOT / "work" / "external_data" / "esm_v3" / "esm_v3_external_validation.csv"
DATA_META = DATA_FILE.with_suffix(".metadata.json")
SUPPLEMENT = core.SUPPLEMENT_DIR
FIGURES = core.FIGURE_DIR

OUT_AUDIT = SUPPLEMENT / "esm_external_validation.md"
OUT_SELECTION = SUPPLEMENT / "esm_external_selection_audit.csv"
OUT_DECOMPOSITION = SUPPLEMENT / "esm_external_decomposition_metrics.csv"
OUT_TERMS = SUPPLEMENT / "esm_external_station_terms.csv"
OUT_REPEATABILITY = SUPPLEMENT / "esm_external_repeatability_metrics.csv"
OUT_REPEATABILITY_PAIRS = SUPPLEMENT / "esm_external_repeatability_pairs.csv"
OUT_FIELD_COMPARISON = SUPPLEMENT / "esm_japan_station_field_comparison.csv"
OUT_TRANSFER_METRICS = SUPPLEMENT / "esm_cross_region_transfer_metrics.csv"
OUT_TRANSFER_PREDICTIONS = SUPPLEMENT / "esm_cross_region_transfer_predictions.csv"
OUT_SPATIAL_METRICS = SUPPLEMENT / "esm_local_spatial_prediction_metrics.csv"
OUT_SPATIAL_PREDICTIONS = SUPPLEMENT / "esm_local_spatial_prediction_predictions.csv"
OUT_FIGURE_PDF = FIGURES / "figure_esm_external_validation.pdf"
OUT_FIGURE_PNG = FIGURES / "figure_esm_external_validation.png"

PERIODS = (1.0, 2.0, 3.0)
PERIOD_SUFFIX = {1.0: "1_000", 2.0: "2_000", 3.0: "3_000"}
MAXIMUM_MW = {"bindi2014_rhyp": 7.6, "kotha2020_rjb": 7.4}
SEED = 20_260_711
REQUESTED_SPATIAL_BLOCKS = 5
MINIMUM_SPATIAL_BLOCK_STATIONS = 30

MODEL_FEATURES = {
    "common_site_hgb": (
        ["log_vs30", "st_elevation"],
        [],
    ),
    "esm_metadata_hgb": (
        ["log_vs30", "st_elevation", "slope_deg"],
        [
            "preferred_estimation_method_vs30_ec8",
            "preferred_ec8_code",
            "reference_site",
            "housing",
            "installation",
        ],
    ),
    "esm_metadata_location_hgb": (
        ["log_vs30", "st_elevation", "slope_deg", "st_longitude", "st_latitude"],
        [
            "preferred_estimation_method_vs30_ec8",
            "preferred_ec8_code",
            "reference_site",
            "housing",
            "installation",
        ],
    ),
}

QUERY_FIELDS = [
    "esm_event_id", "event_time", "ev_nation_code", "ev_latitude", "ev_longitude",
    "ev_depth_km", "fm_type_code", "mw", "ml", "m", "network_code",
    "station_code", "location_code", "channel_code", "sensor_depth_m", "proximity",
    "housing", "installation", "reference_site",
    "st_nation_code", "st_latitude", "st_longitude", "st_elevation",
    "preferred_estimation_method_vs30_ec8", "preferred_ec8_code",
    "preferred_vs30_m_s", "vs30_m_s", "ec8_code", "ec8_code_from_geology",
    "slope_deg", "epi_dist", "hyp_dist", "jb_dist", "rup_dist", "processing_status",
    "quality_class", "processing_type",
]
for component in ("u", "v", "rotd50", "rotd100"):
    QUERY_FIELDS.extend(f"{component}_t{suffix}" for suffix in PERIOD_SUFFIX.values())

QUERY = {
    "min-magnitude": "4.5",
    "max-epicentral-distance": "300",
    "channel": "HN,HG,HL",
    "quality-class": "BEST,GOOD",
    "late-triggered": "N",
    "spectra": "SA",
    "include-fields": ",".join(QUERY_FIELDS),
}
QUERY_URL = "https://esm-db.eu/esmws/flatfile/1/query?" + urllib.parse.urlencode(QUERY)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_data(path: Path, force: bool) -> None:
    if path.exists() and not force:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    request = urllib.request.Request(QUERY_URL, headers={"User-Agent": "japan-station-terms/1.0"})
    with urllib.request.urlopen(request, timeout=240) as response, temporary.open("wb") as output:
        while block := response.read(1024 * 1024):
            output.write(block)
    temporary.replace(path)


def inventory(frame: pd.DataFrame, analysis: str, stage: str) -> dict[str, int | str]:
    return {
        "analysis": analysis,
        "stage": stage,
        "n_records": len(frame),
        "n_events": frame["esm_event_id"].nunique() if "esm_event_id" in frame else 0,
        "n_stations": frame["station_key"].nunique() if "station_key" in frame else 0,
    }


def add_station_key(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    location = frame["location_code"].astype("string").fillna("--")
    frame["station_key"] = (
        frame["network_code"].astype("string").fillna("--")
        + "." + frame["station_code"].astype("string").fillna("--")
        + "." + location
    )
    frame["event_time"] = pd.to_datetime(frame["event_time"], errors="coerce")
    return frame


def deduplicate(frame: pd.DataFrame) -> pd.DataFrame:
    ranked = frame.copy()
    ranked["_quality_rank"] = ranked["quality_class"].map({"BEST": 0, "GOOD": 1}).fillna(9)
    ranked["_processing_rank"] = np.where(
        ranked["processing_type"].astype("string").str.startswith("manual", na=False), 0, 1
    )
    ranked["_channel_rank"] = ranked["channel_code"].map({"HN": 0, "HG": 1, "HL": 2}).fillna(9)
    return (
        ranked.sort_values(["_quality_rank", "_processing_rank", "_channel_rank"])
        .drop_duplicates(["esm_event_id", "station_key"])
        .drop(columns=["_quality_rank", "_processing_rank", "_channel_rank"])
    )


def enforce_support(frame: pd.DataFrame, minimum_event_records: int = 5,
                    minimum_station_records: int = 3) -> pd.DataFrame:
    supported = frame.copy()
    for _ in range(20):
        before = len(supported)
        event_counts = supported["esm_event_id"].value_counts()
        station_counts = supported["station_key"].value_counts()
        supported = supported[
            supported["esm_event_id"].isin(event_counts[event_counts >= minimum_event_records].index)
            & supported["station_key"].isin(station_counts[station_counts >= minimum_station_records].index)
        ]
        if len(supported) == before:
            return supported.copy()
    raise RuntimeError("ESM support filtering did not converge")


def select_records(raw: pd.DataFrame, backbone: str, minimum_mw: float,
                   start_date: str | None) -> tuple[pd.DataFrame, list[dict[str, int | str]]]:
    maximum_mw = MAXIMUM_MW[backbone]
    label = f"{backbone}_mw{minimum_mw:g}-{maximum_mw:g}_{start_date or 'all'}"
    frame = add_station_key(raw)
    audit = [inventory(frame, label, "raw_query")]

    common = (
        frame["mw"].between(minimum_mw, maximum_mw)
        & frame["ev_depth_km"].between(0, 35)
        & frame["preferred_vs30_m_s"].gt(0)
        & frame["st_latitude"].between(25, 72)
        & frame["st_longitude"].between(-25, 45)
        & frame["ev_latitude"].between(25, 72)
        & frame["ev_longitude"].between(-25, 45)
        & pd.to_numeric(frame["sensor_depth_m"], errors="coerce").fillna(0).le(0.5)
        & ~frame["proximity"].isin(["Inside structure", "Close to structure"])
        & frame["quality_class"].isin(["BEST", "GOOD"])
    )
    if start_date:
        common &= frame["event_time"].ge(start_date)
    if backbone == "bindi2014_rhyp":
        common &= frame["hyp_dist"].between(1, 300)
        for suffix in PERIOD_SUFFIX.values():
            common &= frame[f"u_t{suffix}"].gt(0) & frame[f"v_t{suffix}"].gt(0)
    elif backbone == "kotha2020_rjb":
        common &= frame["jb_dist"].between(0, 300)
        for suffix in PERIOD_SUFFIX.values():
            common &= frame[f"rotd50_t{suffix}"].gt(0)
    else:
        raise ValueError(backbone)

    selected = frame[common].copy()
    audit.append(inventory(selected, label, "physical_and_response_filters"))
    selected = deduplicate(selected)
    audit.append(inventory(selected, label, "event_station_deduplication"))
    selected = enforce_support(selected)
    audit.append(inventory(selected, label, "iterative_event_station_support"))

    station_ids = {key: index + 1 for index, key in enumerate(sorted(selected["station_key"].unique()))}
    event_ids = {key: index + 1 for index, key in enumerate(sorted(selected["esm_event_id"].unique()))}
    selected["siteid2"] = selected["station_key"].map(station_ids).astype(int)
    selected["eq_source_id"] = selected["esm_event_id"].map(event_ids).astype(int)
    return selected.reset_index(drop=True), audit


def residual_frames(records: pd.DataFrame, backbone: str) -> dict[float, pd.DataFrame]:
    n_records = len(records)
    if backbone == "bindi2014_rhyp":
        model = BindiEtAl2014RhypEC8NoSOF()
        context = np.rec.fromarrays(
            [
                records["preferred_vs30_m_s"].to_numpy(float),
                records["mw"].to_numpy(float),
                np.zeros(n_records),
                records["hyp_dist"].to_numpy(float),
            ],
            names=["vs30", "mag", "rake", "rhypo"],
        )
    else:
        model = KothaEtAl2020Site()
        context = np.rec.fromarrays(
            [
                records["preferred_vs30_m_s"].to_numpy(float),
                records["mw"].to_numpy(float),
                records["ev_depth_km"].to_numpy(float),
                records["jb_dist"].to_numpy(float),
            ],
            names=["vs30", "mag", "hypo_depth", "rjb"],
        )

    mean = np.zeros((len(PERIODS), n_records))
    zeros = np.zeros_like(mean)
    model.compute(context, [SA(period) for period in PERIODS], mean, zeros.copy(), zeros.copy(), zeros.copy())

    output: dict[float, pd.DataFrame] = {}
    metadata = [
        "siteid2", "eq_source_id", "station_key", "esm_event_id", "event_time",
        "network_code", "station_code", "location_code", "st_nation_code", "st_latitude",
        "st_longitude", "st_elevation", "preferred_vs30_m_s", "slope_deg",
        "preferred_estimation_method_vs30_ec8", "preferred_ec8_code", "reference_site",
        "housing", "installation",
    ]
    for index, period in enumerate(PERIODS):
        suffix = PERIOD_SUFFIX[period]
        block = records[metadata].copy()
        if backbone == "bindi2014_rhyp":
            observed_g = np.sqrt(records[f"u_t{suffix}"] * records[f"v_t{suffix}"]) / core.G_IN_CM_S2
            block["response_component"] = "geometric_mean_as_recorded"
        else:
            observed_g = records[f"rotd50_t{suffix}"] / core.G_IN_CM_S2
            block["response_component"] = "RotD50"
        block["residual_site"] = np.log10(observed_g.to_numpy(float)) - mean[index] / math.log(10)
        output[period] = block
    return output


def full_station_terms(frames: dict[float, pd.DataFrame], backbone: str,
                       minimum_mw: float, start_date: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    term_blocks: list[pd.DataFrame] = []
    metric_rows: list[dict[str, float | int | str]] = []
    for period, frame in frames.items():
        terms, components = component_station_terms(frame)
        metadata = frame.sort_values("station_key").drop_duplicates("siteid2")[
            ["siteid2", "station_key", "network_code", "station_code", "location_code",
             "st_nation_code", "st_latitude", "st_longitude", "st_elevation",
             "preferred_vs30_m_s", "slope_deg", "preferred_estimation_method_vs30_ec8",
             "preferred_ec8_code", "reference_site", "housing", "installation",
             "response_component"]
        ]
        terms = terms.merge(metadata, on="siteid2", validate="one_to_one")
        terms["period_s"] = period
        terms["backbone"] = backbone
        terms["minimum_mw"] = minimum_mw
        terms["start_date"] = start_date or "all"
        term_blocks.append(terms)
        metric_rows.append(
            {
                "backbone": backbone,
                "minimum_mw": minimum_mw,
                "start_date": start_date or "all",
                "period_s": period,
                "n_records": len(frame),
                "n_events": frame["eq_source_id"].nunique(),
                "n_stations": frame["siteid2"].nunique(),
                "n_graph_components": len(components),
                "largest_component_records": int(components["component_records"].max()),
                "solver_max_relative_normal_residual": float(
                    components["solver_relative_normal_residual"].max()
                ),
            }
        )
    return pd.concat(term_blocks, ignore_index=True), pd.DataFrame(metric_rows)


def center_by_group(frame: pd.DataFrame) -> pd.DataFrame:
    centered = frame.copy()
    centered["train_centered_log10"] = np.nan
    centered["test_centered_log10"] = np.nan
    for _, block in centered.groupby("component_pair"):
        index = block.index
        weights = block["n_records_test"].to_numpy(float)
        centered.loc[index, "train_centered_log10"] = (
            block["station_effect_log10_train"] - np.average(
                block["station_effect_log10_train"], weights=weights
            )
        )
        centered.loc[index, "test_centered_log10"] = (
            block["station_effect_log10_test"] - np.average(
                block["station_effect_log10_test"], weights=weights
            )
        )
    return centered


def repeatability(frames: dict[float, pd.DataFrame], backbone: str, minimum_mw: float,
                  start_date: str | None, minimum_train_records: int,
                  minimum_test_records: int, minimum_component_pair_stations: int,
                  seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, float | int | str]] = []
    pair_blocks: list[pd.DataFrame] = []
    for period, frame in frames.items():
        events = np.array(sorted(frame["eq_source_id"].unique()))
        rng = np.random.default_rng(seed)
        rng.shuffle(events)
        fold_map = {int(event): index % 5 for index, event in enumerate(events)}
        working = frame.copy()
        working["fold"] = working["eq_source_id"].map(fold_map).astype(int)
        for fold in range(5):
            train, train_components = component_station_terms(working[working["fold"].ne(fold)])
            test, test_components = component_station_terms(working[working["fold"].eq(fold)])
            paired = train[train["n_records"].ge(minimum_train_records)].merge(
                test[test["n_records"].ge(minimum_test_records)],
                on="siteid2", suffixes=("_train", "_test"), validate="one_to_one",
            )
            paired["component_pair"] = (
                paired["graph_component_train"].astype(str)
                + "|" + paired["graph_component_test"].astype(str)
            )
            pair_size = paired.groupby("component_pair")["siteid2"].transform("size")
            paired = center_by_group(paired[pair_size.ge(minimum_component_pair_stations)].copy())
            if len(paired) < 3:
                continue
            weights = paired["n_records_test"].to_numpy(float)
            train_values = paired["train_centered_log10"].to_numpy(float)
            test_values = paired["test_centered_log10"].to_numpy(float)
            baseline_rmse = core.weighted_rmse(test_values, np.zeros_like(test_values), weights)
            prediction_rmse = core.weighted_rmse(test_values, train_values, weights)
            metric_rows.append(
                {
                    "backbone": backbone,
                    "minimum_mw": minimum_mw,
                    "start_date": start_date or "all",
                    "period_s": period,
                    "fold": fold,
                    "n_train_events": working.loc[working["fold"].ne(fold), "eq_source_id"].nunique(),
                    "n_test_events": working.loc[working["fold"].eq(fold), "eq_source_id"].nunique(),
                    "n_paired_stations": len(paired),
                    "n_component_pairs": paired["component_pair"].nunique(),
                    "n_train_graph_components": len(train_components),
                    "n_test_graph_components": len(test_components),
                    "pearson": float(np.corrcoef(train_values, test_values)[0, 1]),
                    "zero_baseline_rmse_log10": baseline_rmse,
                    "train_term_rmse_log10": prediction_rmse,
                    "rmse_gain_pct": 100.0 * (1.0 - prediction_rmse / baseline_rmse),
                }
            )
            paired["backbone"] = backbone
            paired["minimum_mw"] = minimum_mw
            paired["start_date"] = start_date or "all"
            paired["period_s"] = period
            paired["fold"] = fold
            pair_blocks.append(paired)
    return pd.DataFrame(metric_rows), pd.concat(pair_blocks, ignore_index=True)


def bootstrap_repeatability(pairs: pd.DataFrame, replicates: int, seed: int) -> tuple[float, float, float, float]:
    stations = np.array(sorted(pairs["siteid2"].unique()))
    station_codes = pd.Categorical(pairs["siteid2"], categories=stations).codes
    rng = np.random.default_rng(seed)
    correlations = np.empty(replicates)
    gains = np.empty(replicates)
    for replicate in range(replicates):
        counts = np.bincount(rng.integers(0, len(stations), len(stations)), minlength=len(stations))
        sampled = counts[station_codes]
        fold_correlations = []
        fold_gains = []
        for _, fold in pairs[sampled.gt(0) if isinstance(sampled, pd.Series) else sampled > 0].groupby("fold"):
            index = fold.index.to_numpy()
            multiplicity = sampled[index].astype(float)
            x = fold["train_centered_log10"].to_numpy(float)
            y = fold["test_centered_log10"].to_numpy(float)
            record_weights = fold["n_records_test"].to_numpy(float) * multiplicity
            mean_x = np.average(x, weights=multiplicity)
            mean_y = np.average(y, weights=multiplicity)
            covariance = np.average((x - mean_x) * (y - mean_y), weights=multiplicity)
            denominator = np.sqrt(
                np.average((x - mean_x) ** 2, weights=multiplicity)
                * np.average((y - mean_y) ** 2, weights=multiplicity)
            )
            if denominator > np.finfo(float).eps:
                fold_correlations.append(covariance / denominator)
            baseline = core.weighted_rmse(y, np.zeros_like(y), record_weights)
            prediction = core.weighted_rmse(y, x, record_weights)
            if baseline > np.finfo(float).eps:
                fold_gains.append(100.0 * (1.0 - prediction / baseline))
        correlations[replicate] = np.mean(fold_correlations) if fold_correlations else np.nan
        gains[replicate] = np.mean(fold_gains) if fold_gains else np.nan
    corr_ci = np.nanquantile(correlations, [0.025, 0.975])
    gain_ci = np.nanquantile(gains, [0.025, 0.975])
    return float(corr_ci[0]), float(corr_ci[1]), float(gain_ci[0]), float(gain_ci[1])


def aggregate_repeatability(metrics: pd.DataFrame, pairs: pd.DataFrame,
                            replicates: int, seed: int) -> pd.DataFrame:
    rows = []
    keys = ["backbone", "minimum_mw", "start_date", "period_s"]
    for key, block in metrics.groupby(keys, sort=False):
        matching = pairs.copy()
        for column, value in zip(keys, key, strict=True):
            matching = matching[matching[column].eq(value)]
        corr_low, corr_high, gain_low, gain_high = bootstrap_repeatability(
            matching.reset_index(drop=True), replicates, seed + int(round(100 * key[-1]))
        )
        rows.append(
            {
                **dict(zip(keys, key, strict=True)),
                "scope": "fold_mean",
                "n_folds": len(block),
                "mean_paired_stations": block["n_paired_stations"].mean(),
                "pearson": block["pearson"].mean(),
                "pearson_ci_low": corr_low,
                "pearson_ci_high": corr_high,
                "rmse_gain_pct": block["rmse_gain_pct"].mean(),
                "rmse_gain_pct_ci_low": gain_low,
                "rmse_gain_pct_ci_high": gain_high,
            }
        )
    return pd.concat([metrics.assign(scope="fold"), pd.DataFrame(rows)], ignore_index=True)


def weighted_summary(values: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    order = np.argsort(values)
    ordered = values[order]
    cumulative = np.cumsum(weights[order]) / np.sum(weights)
    quantiles = np.interp([0.05, 0.5, 0.95], cumulative, ordered)
    centered = values - np.average(values, weights=weights)
    return {
        "weighted_std_log10": float(np.sqrt(np.average(centered**2, weights=weights))),
        "multiplier_p05": float(10.0 ** quantiles[0]),
        "multiplier_median": float(10.0 ** quantiles[1]),
        "multiplier_p95": float(10.0 ** quantiles[2]),
    }


def field_comparison(esm_terms: pd.DataFrame, minimum_records: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    esm_primary = esm_terms[
        esm_terms["backbone"].eq("bindi2014_rhyp")
        & esm_terms["minimum_mw"].eq(4.5)
        & esm_terms["start_date"].eq("2015-01-01")
        & esm_terms["n_records"].ge(minimum_records)
    ].copy()
    for period, block in esm_primary.groupby("period_s"):
        summary = weighted_summary(
            block["station_effect_log10"].to_numpy(float), block["n_records"].to_numpy(float)
        )
        rows.append({"region": "ESM Europe", "period_s": period, "n_stations": len(block), **summary})

    component = pd.read_csv(SUPPLEMENT / "jshis_response_component_station_terms.csv")
    metadata = pd.read_csv(SUPPLEMENT / "jshis_event_adjusted_station_terms.csv")
    metadata = metadata[metadata["model"].eq("mf2013_site")][
        ["siteid2", "period_s", "avs30", "elevation", "network_label"]
    ]
    japan = component.merge(metadata, on=["siteid2", "period_s"], validate="one_to_one")
    japan = japan[
        japan["period_s"].isin(PERIODS)
        & japan["n_records_rotd50"].ge(minimum_records)
    ].copy()
    for period, block in japan.groupby("period_s"):
        summary = weighted_summary(
            block["station_effect_rotd50_log10"].to_numpy(float),
            block["n_records_rotd50"].to_numpy(float),
        )
        rows.append({"region": "J-SHIS Japan", "period_s": period, "n_stations": len(block), **summary})
    return pd.DataFrame(rows), japan


def bootstrap_transfer(observed: np.ndarray, predicted: np.ndarray, weights: np.ndarray,
                       replicates: int, seed: int) -> tuple[float, float, float, float]:
    rng = np.random.default_rng(seed)
    correlations = np.empty(replicates)
    gains = np.empty(replicates)
    for index in range(replicates):
        sample = rng.integers(0, len(observed), len(observed))
        obs = observed[sample]
        pred = predicted[sample]
        sample_weights = weights[sample]
        correlations[index] = np.corrcoef(obs, pred)[0, 1]
        baseline = core.weighted_rmse(obs, np.zeros_like(obs), sample_weights)
        model = core.weighted_rmse(obs, pred, sample_weights)
        gains[index] = 100.0 * (1.0 - model / baseline)
    corr_ci = np.quantile(correlations, [0.025, 0.975])
    gain_ci = np.quantile(gains, [0.025, 0.975])
    return float(corr_ci[0]), float(corr_ci[1]), float(gain_ci[0]), float(gain_ci[1])


def cross_region_transfer(esm_terms: pd.DataFrame, japan: pd.DataFrame,
                          replicates: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    prediction_blocks = []
    target = esm_terms[
        esm_terms["backbone"].eq("bindi2014_rhyp")
        & esm_terms["minimum_mw"].eq(4.5)
        & esm_terms["start_date"].eq("2015-01-01")
        & esm_terms["n_records"].ge(15)
    ].copy()
    source = japan[
        japan["network_label"].eq("K-NET")
        & japan["n_records_rotd50"].ge(15)
        & japan["avs30"].gt(0)
        & japan["elevation"].notna()
    ].copy()
    source["log_vs30"] = np.log10(source["avs30"])
    target["log_vs30"] = np.log10(target["preferred_vs30_m_s"])
    for period in PERIODS:
        train = source[source["period_s"].eq(period)].copy()
        test = target[target["period_s"].eq(period)].copy()
        model = HistGradientBoostingRegressor(
            loss="squared_error", learning_rate=0.04, max_iter=400, max_leaf_nodes=15,
            min_samples_leaf=35, l2_regularization=0.1, random_state=seed,
        )
        train_x = train[["log_vs30", "elevation"]].to_numpy(float)
        train_y = train["station_effect_rotd50_log10"].to_numpy(float)
        train_weights = train["n_records_rotd50"].to_numpy(float)
        model.fit(train_x, train_y, sample_weight=train_weights)
        train_prediction = model.predict(train_x)
        prediction = model.predict(test[["log_vs30", "st_elevation"]].to_numpy(float))
        prediction -= np.average(train_prediction, weights=train_weights)

        observed = test["station_effect_log10"].to_numpy(float).copy()
        test_weights = test["n_records"].to_numpy(float)
        # The target field and prediction are separately centred to remove the
        # arbitrary additive offset of the station-event decomposition.
        observed -= np.average(observed, weights=test_weights)
        prediction -= np.average(prediction, weights=test_weights)
        correlation = float(np.corrcoef(observed, prediction)[0, 1])
        baseline = core.weighted_rmse(observed, np.zeros_like(observed), test_weights)
        prediction_rmse = core.weighted_rmse(observed, prediction, test_weights)
        gain = 100.0 * (1.0 - prediction_rmse / baseline)
        corr_low, corr_high, gain_low, gain_high = bootstrap_transfer(
            observed, prediction, test_weights, replicates, seed + int(100 * period)
        )
        metric_rows.append(
            {
                "period_s": period,
                "source_region": "Japan K-NET",
                "target_region": "ESM Europe post-2015",
                "features": "log10(Vs30), elevation",
                "prediction_centering": "target_prediction_weighted_zero_without_labels",
                "n_source_stations": len(train),
                "n_target_stations": len(test),
                "pearson": correlation,
                "pearson_ci_low": corr_low,
                "pearson_ci_high": corr_high,
                "rmse_gain_pct": gain,
                "rmse_gain_pct_ci_low": gain_low,
                "rmse_gain_pct_ci_high": gain_high,
            }
        )
        block = test[
            ["station_key", "st_nation_code", "st_latitude", "st_longitude",
             "preferred_vs30_m_s", "st_elevation", "n_records"]
        ].copy()
        block["period_s"] = period
        block["observed_station_term_log10"] = observed
        block["frozen_prediction_log10"] = prediction
        block["prediction_error_log10"] = observed - prediction
        prediction_blocks.append(block)
    return pd.DataFrame(metric_rows), pd.concat(prediction_blocks, ignore_index=True)


def build_local_model(model_name: str, seed: int, training_frame: pd.DataFrame) -> Pipeline:
    numeric_features, categorical_features = MODEL_FEATURES[model_name]
    categorical_features = [
        column for column in categorical_features if training_frame[column].notna().any()
    ]
    transformers = [
        (
            "numeric",
            SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
            numeric_features,
        )
    ]
    if categorical_features:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                    ]
                ),
                categorical_features,
            )
        )
    estimator = HistGradientBoostingRegressor(
        loss="squared_error", learning_rate=0.04, max_iter=400, max_leaf_nodes=15,
        min_samples_leaf=35, l2_regularization=0.1, random_state=seed,
    )
    return Pipeline(
        [("preprocess", ColumnTransformer(transformers)), ("model", estimator)]
    )


def choose_spatial_folds(frame: pd.DataFrame, seed: int) -> tuple[np.ndarray, int]:
    coordinates = StandardScaler().fit_transform(
        frame[["st_longitude", "st_latitude"]].to_numpy(float)
    )
    maximum = min(REQUESTED_SPATIAL_BLOCKS, len(frame))
    for n_splits in range(maximum, 1, -1):
        folds = KMeans(n_clusters=n_splits, random_state=seed, n_init=30).fit_predict(coordinates)
        if np.bincount(folds, minlength=n_splits).min() >= MINIMUM_SPATIAL_BLOCK_STATIONS:
            return folds, n_splits
    raise ValueError("ESM station geometry cannot support two spatial blocks")


def prediction_statistics(observed: np.ndarray, predicted: np.ndarray,
                          weights: np.ndarray) -> dict[str, float]:
    baseline = core.weighted_rmse(observed, np.zeros_like(observed), weights)
    model_rmse = core.weighted_rmse(observed, predicted, weights)
    return {
        "pearson": float(np.corrcoef(observed, predicted)[0, 1]),
        "zero_baseline_rmse_log10": baseline,
        "prediction_rmse_log10": model_rmse,
        "rmse_gain_pct": 100.0 * (1.0 - model_rmse / baseline),
    }


def multiplicity_correlation(observed: np.ndarray, predicted: np.ndarray,
                             multiplicity: np.ndarray) -> float:
    observed_mean = np.average(observed, weights=multiplicity)
    predicted_mean = np.average(predicted, weights=multiplicity)
    covariance = np.average(
        (observed - observed_mean) * (predicted - predicted_mean), weights=multiplicity
    )
    denominator = np.sqrt(
        np.average((observed - observed_mean) ** 2, weights=multiplicity)
        * np.average((predicted - predicted_mean) ** 2, weights=multiplicity)
    )
    return float(covariance / denominator) if denominator > np.finfo(float).eps else np.nan


def bootstrap_prediction_statistics(predictions: pd.DataFrame, replicates: int,
                                    seed: int, group_column: str | None) -> tuple[float, ...]:
    stations = np.array(sorted(predictions["siteid2"].unique()))
    station_codes = pd.Categorical(predictions["siteid2"], categories=stations).codes
    groups = [None] if group_column is None else sorted(predictions[group_column].unique())
    rng = np.random.default_rng(seed)
    correlations = np.full(replicates, np.nan)
    gains = np.full(replicates, np.nan)
    for replicate in range(replicates):
        counts = np.bincount(
            rng.integers(0, len(stations), len(stations)), minlength=len(stations)
        )
        multiplicity = counts[station_codes].astype(float)
        group_correlations = []
        group_gains = []
        for group in groups:
            selected = multiplicity > 0
            if group is not None:
                selected &= predictions[group_column].to_numpy() == group
            block = predictions.loc[selected]
            sampled = multiplicity[selected]
            observed = block["observed_station_term_log10"].to_numpy(float)
            predicted = block["oof_prediction_log10"].to_numpy(float)
            weights = block["test_record_weight"].to_numpy(float) * sampled
            correlation = multiplicity_correlation(observed, predicted, sampled)
            if np.isfinite(correlation):
                group_correlations.append(correlation)
            baseline = core.weighted_rmse(observed, np.zeros_like(observed), weights)
            model_rmse = core.weighted_rmse(observed, predicted, weights)
            if baseline > np.finfo(float).eps:
                group_gains.append(100.0 * (1.0 - model_rmse / baseline))
        if group_correlations:
            correlations[replicate] = np.mean(group_correlations)
        if group_gains:
            gains[replicate] = np.mean(group_gains)
    corr_ci = np.nanquantile(correlations, [0.025, 0.975])
    gain_ci = np.nanquantile(gains, [0.025, 0.975])
    return float(corr_ci[0]), float(corr_ci[1]), float(gain_ci[0]), float(gain_ci[1])


def primary_spatial_frame(terms: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, int]:
    frame = terms[
        terms["backbone"].eq("bindi2014_rhyp")
        & terms["minimum_mw"].eq(4.5)
        & terms["start_date"].eq("2015-01-01")
        & terms["n_records"].ge(15)
    ].copy()
    reference = frame[frame["period_s"].eq(PERIODS[0])]
    largest_component = reference.groupby("graph_component")["n_records"].sum().idxmax()
    frame = frame[frame["graph_component"].eq(largest_component)].copy()
    frame["log_vs30"] = np.log10(frame["preferred_vs30_m_s"])
    reference = frame[frame["period_s"].eq(PERIODS[0])].sort_values("siteid2")
    folds, n_splits = choose_spatial_folds(reference, seed)
    fold_map = dict(zip(reference["siteid2"], folds, strict=True))
    frame["spatial_fold"] = frame["siteid2"].map(fold_map).astype(int)
    if frame.groupby("period_s")["siteid2"].nunique().nunique() != 1:
        raise ValueError("ESM primary station set differs between periods")
    return frame, n_splits


def full_field_spatial_prediction(frame: pd.DataFrame, n_splits: int, replicates: int,
                                  seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    prediction_blocks = []
    for period in PERIODS:
        period_frame = frame[frame["period_s"].eq(period)].copy().reset_index(drop=True)
        observed = period_frame["station_effect_log10"].to_numpy(float)
        weights = period_frame["n_records"].to_numpy(float)
        folds = period_frame["spatial_fold"].to_numpy(int)
        for model_name in MODEL_FEATURES:
            raw_prediction = np.full(len(period_frame), np.nan)
            for fold in range(n_splits):
                train = folds != fold
                test = folds == fold
                train_target = observed[train] - np.average(observed[train], weights=weights[train])
                model = build_local_model(model_name, seed + fold, period_frame.loc[train])
                model.fit(
                    period_frame.loc[train], train_target,
                    model__sample_weight=weights[train],
                )
                raw_prediction[test] = model.predict(period_frame.loc[test])
            prediction = raw_prediction - np.average(raw_prediction, weights=weights)
            block = period_frame[
                ["siteid2", "station_key", "st_nation_code", "st_latitude", "st_longitude",
                 "preferred_vs30_m_s", "st_elevation", "slope_deg", "n_records",
                 "spatial_fold"]
            ].copy()
            block["validation_design"] = "full_field_spatial_blocks"
            block["event_fold"] = -1
            block["model"] = model_name
            block["period_s"] = period
            block["observed_station_term_log10"] = observed
            block["raw_oof_prediction_log10"] = raw_prediction
            block["oof_prediction_log10"] = prediction
            block["test_record_weight"] = weights
            prediction_blocks.append(block)
            for fold in range(n_splits):
                selected = folds == fold
                metric_rows.append(
                    {
                        "validation_design": "full_field_spatial_blocks",
                        "period_s": period,
                        "model": model_name,
                        "scope": "spatial_fold",
                        "event_fold": -1,
                        "spatial_fold": fold,
                        "n_spatial_blocks": n_splits,
                        "n_stations": int(selected.sum()),
                        "n_records_weight": int(weights[selected].sum()),
                        **prediction_statistics(observed[selected], prediction[selected], weights[selected]),
                    }
                )
            statistics = prediction_statistics(observed, prediction, weights)
            corr_low, corr_high, gain_low, gain_high = bootstrap_prediction_statistics(
                block, replicates, seed + int(100 * period), None
            )
            metric_rows.append(
                {
                    "validation_design": "full_field_spatial_blocks",
                    "period_s": period,
                    "model": model_name,
                    "scope": "overall",
                    "event_fold": -1,
                    "spatial_fold": -1,
                    "n_spatial_blocks": n_splits,
                    "n_stations": len(period_frame),
                    "n_records_weight": int(weights.sum()),
                    **statistics,
                    "pearson_ci_low": corr_low,
                    "pearson_ci_high": corr_high,
                    "rmse_gain_pct_ci_low": gain_low,
                    "rmse_gain_pct_ci_high": gain_high,
                }
            )
    return pd.DataFrame(metric_rows), pd.concat(prediction_blocks, ignore_index=True)


def held_event_spatial_prediction(frame: pd.DataFrame, pairs: pd.DataFrame, n_splits: int,
                                  replicates: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary_pairs = pairs[
        pairs["backbone"].eq("bindi2014_rhyp")
        & pairs["minimum_mw"].eq(4.5)
        & pairs["start_date"].eq("2015-01-01")
    ].copy()
    metadata_columns = [
        "siteid2", "period_s", "station_key", "st_nation_code", "st_latitude",
        "st_longitude", "st_elevation", "preferred_vs30_m_s", "slope_deg",
        "preferred_estimation_method_vs30_ec8", "preferred_ec8_code", "reference_site",
        "housing", "installation", "spatial_fold", "log_vs30",
    ]
    primary_pairs = primary_pairs.merge(
        frame[metadata_columns], on=["siteid2", "period_s"], validate="many_to_one"
    )
    metric_rows = []
    prediction_blocks = []
    for (period, event_fold), period_frame in primary_pairs.groupby(
        ["period_s", "fold"], sort=True
    ):
        period_frame = period_frame.copy().reset_index(drop=True)
        observed = period_frame["test_centered_log10"].to_numpy(float).copy()
        test_weights = period_frame["n_records_test"].to_numpy(float)
        observed -= np.average(observed, weights=test_weights)
        training_field = period_frame["train_centered_log10"].to_numpy(float)
        train_weights = period_frame["n_records_train"].to_numpy(float)
        spatial_folds = period_frame["spatial_fold"].to_numpy(int)
        for model_name in MODEL_FEATURES:
            raw_prediction = np.full(len(period_frame), np.nan)
            for spatial_fold in range(n_splits):
                train = spatial_folds != spatial_fold
                test = spatial_folds == spatial_fold
                if not test.any():
                    continue
                train_target = training_field[train] - np.average(
                    training_field[train], weights=train_weights[train]
                )
                model = build_local_model(
                    model_name, seed + 10 * int(event_fold) + spatial_fold,
                    period_frame.loc[train],
                )
                model.fit(
                    period_frame.loc[train], train_target,
                    model__sample_weight=train_weights[train],
                )
                raw_prediction[test] = model.predict(period_frame.loc[test])
            if np.isnan(raw_prediction).any():
                raise ValueError("A held-event ESM station was not assigned an out-of-fold prediction")
            prediction = raw_prediction - np.average(raw_prediction, weights=test_weights)
            statistics = prediction_statistics(observed, prediction, test_weights)
            metric_rows.append(
                {
                    "validation_design": "held_event_field_spatial_blocks",
                    "period_s": period,
                    "model": model_name,
                    "scope": "event_fold",
                    "event_fold": int(event_fold),
                    "spatial_fold": -1,
                    "n_spatial_blocks": n_splits,
                    "n_stations": len(period_frame),
                    "n_records_weight": int(test_weights.sum()),
                    **statistics,
                }
            )
            block = period_frame[
                ["siteid2", "station_key", "st_nation_code", "st_latitude", "st_longitude",
                 "preferred_vs30_m_s", "st_elevation", "slope_deg", "spatial_fold"]
            ].copy()
            block["validation_design"] = "held_event_field_spatial_blocks"
            block["event_fold"] = int(event_fold)
            block["model"] = model_name
            block["period_s"] = period
            block["observed_station_term_log10"] = observed
            block["raw_oof_prediction_log10"] = raw_prediction
            block["oof_prediction_log10"] = prediction
            block["test_record_weight"] = test_weights
            prediction_blocks.append(block)

    predictions = pd.concat(prediction_blocks, ignore_index=True)
    fold_metrics = pd.DataFrame(metric_rows)
    aggregate_rows = []
    for (period, model_name), block in fold_metrics.groupby(["period_s", "model"], sort=True):
        selected_predictions = predictions[
            predictions["period_s"].eq(period) & predictions["model"].eq(model_name)
        ]
        corr_low, corr_high, gain_low, gain_high = bootstrap_prediction_statistics(
            selected_predictions, replicates, seed + int(100 * period), "event_fold"
        )
        aggregate_rows.append(
            {
                "validation_design": "held_event_field_spatial_blocks",
                "period_s": period,
                "model": model_name,
                "scope": "event_fold_mean",
                "event_fold": -1,
                "spatial_fold": -1,
                "n_spatial_blocks": n_splits,
                "n_stations": int(selected_predictions["siteid2"].nunique()),
                "n_records_weight": int(selected_predictions["test_record_weight"].sum()),
                "pearson": block["pearson"].mean(),
                "zero_baseline_rmse_log10": block["zero_baseline_rmse_log10"].mean(),
                "prediction_rmse_log10": block["prediction_rmse_log10"].mean(),
                "rmse_gain_pct": block["rmse_gain_pct"].mean(),
                "pearson_ci_low": corr_low,
                "pearson_ci_high": corr_high,
                "rmse_gain_pct_ci_low": gain_low,
                "rmse_gain_pct_ci_high": gain_high,
            }
        )
    return pd.concat([fold_metrics, pd.DataFrame(aggregate_rows)], ignore_index=True), predictions


def local_spatial_prediction(terms: pd.DataFrame, pairs: pd.DataFrame, replicates: int,
                             seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame, n_splits = primary_spatial_frame(terms, seed)
    full_metrics, full_predictions = full_field_spatial_prediction(
        frame, n_splits, replicates, seed
    )
    strict_metrics, strict_predictions = held_event_spatial_prediction(
        frame, pairs, n_splits, replicates, seed
    )
    return (
        pd.concat([full_metrics, strict_metrics], ignore_index=True),
        pd.concat([full_predictions, strict_predictions], ignore_index=True),
    )


def plot_results(terms: pd.DataFrame, repeatability_metrics: pd.DataFrame,
                 repeatability_pairs: pd.DataFrame, comparison: pd.DataFrame) -> None:
    primary_terms = terms[
        terms["backbone"].eq("bindi2014_rhyp")
        & terms["minimum_mw"].eq(4.5)
        & terms["start_date"].eq("2015-01-01")
        & terms["period_s"].eq(3.0)
        & terms["n_records"].ge(15)
    ].copy()
    aggregate = repeatability_metrics[repeatability_metrics["scope"].eq("fold_mean")]
    primary_pairs = repeatability_pairs[
        repeatability_pairs["backbone"].eq("bindi2014_rhyp")
        & repeatability_pairs["minimum_mw"].eq(4.5)
        & repeatability_pairs["start_date"].eq("2015-01-01")
        & repeatability_pairs["period_s"].eq(3.0)
    ].copy()

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.0), constrained_layout=True)
    ax = axes[0, 0]
    points = ax.scatter(
        primary_terms["st_longitude"], primary_terms["st_latitude"],
        c=primary_terms["station_effect_log10"], s=np.clip(primary_terms["n_records"], 10, 80),
        cmap="RdBu_r", vmin=-0.45, vmax=0.45, alpha=0.8, linewidths=0,
    )
    ax.set(xlabel="Longitude (degrees E)", ylabel="Latitude (degrees N)",
           title="a  Post-2015 ESM station terms at 3.0 s")
    fig.colorbar(points, ax=ax, label="Station term (log$_{10}$)")

    ax = axes[0, 1]
    styles = {
        ("bindi2014_rhyp", 4.5, "2015-01-01"): ("o-", "#0072B2", "Bindi, post-2015"),
        ("bindi2014_rhyp", 5.0, "2015-01-01"): ("s--", "#009E73", "Bindi, $M_w \\geq 5$"),
        ("kotha2020_rjb", 4.5, "all"): ("^-.", "#D55E00", "Kotha RotD50"),
    }
    for key, (style, color, label) in styles.items():
        block = aggregate[
            aggregate["backbone"].eq(key[0])
            & aggregate["minimum_mw"].eq(key[1])
            & aggregate["start_date"].eq(key[2])
        ].sort_values("period_s")
        if block.empty:
            continue
        ax.errorbar(
            block["period_s"], block["pearson"],
            yerr=np.vstack([block["pearson"] - block["pearson_ci_low"],
                            block["pearson_ci_high"] - block["pearson"]]),
            fmt=style, color=color, capsize=3, label=label,
        )
    ax.set(xlabel="Period (s)", ylabel="Held-event station correlation",
           title="b  Independent-event repeatability", xticks=list(PERIODS), ylim=(0, 1.02))
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.scatter(
        primary_pairs["train_centered_log10"], primary_pairs["test_centered_log10"],
        s=np.clip(primary_pairs["n_records_test"], 8, 55), alpha=0.32,
        color="#2878B5", linewidths=0,
    )
    limits = np.quantile(np.r_[primary_pairs["train_centered_log10"],
                                primary_pairs["test_centered_log10"]], [0.005, 0.995])
    ax.plot(limits, limits, "--", color="#444444", lw=1)
    ax.text(0.96, 0.04, "Dashed: equality", transform=ax.transAxes, ha="right", color="#444444")
    row = aggregate[
        aggregate["backbone"].eq("bindi2014_rhyp")
        & aggregate["minimum_mw"].eq(4.5)
        & aggregate["start_date"].eq("2015-01-01")
        & aggregate["period_s"].eq(3.0)
    ].iloc[0]
    ax.text(0.04, 0.96, f"Mean $r$ = {row['pearson']:.3f}\nRMSE gain = {row['rmse_gain_pct']:.1f}%",
            transform=ax.transAxes, va="top")
    ax.set(xlim=limits, ylim=limits, xlabel="Training-event station term (log$_{10}$)",
           ylabel="Held-event station term (log$_{10}$)",
           title="c  Component-centred 3.0 s pairs")

    ax = axes[1, 1]
    colors = {"ESM Europe": "#D55E00", "J-SHIS Japan": "#0072B2"}
    for region, block in comparison.groupby("region"):
        block = block.sort_values("period_s")
        ax.plot(block["period_s"], block["weighted_std_log10"], "o-",
                color=colors[region], label=region)
    ax.set(xlabel="Period (s)", ylabel="Weighted station-field SD (log$_{10}$)",
           title="d  Long-period field amplitude", xticks=list(PERIODS))
    ax.legend(frameon=False)

    for axis in axes.flat:
        axis.grid(alpha=0.2, linewidth=0.6)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audit(raw: pd.DataFrame, metadata: dict, selection: pd.DataFrame,
                terms: pd.DataFrame, repeatability_metrics: pd.DataFrame,
                transfer_metrics: pd.DataFrame, spatial_metrics: pd.DataFrame) -> None:
    aggregate = repeatability_metrics[repeatability_metrics["scope"].eq("fold_mean")]
    primary = aggregate[
        aggregate["backbone"].eq("bindi2014_rhyp")
        & aggregate["minimum_mw"].eq(4.5)
        & aggregate["start_date"].eq("2015-01-01")
    ].sort_values("period_s")
    sensitivity = aggregate[
        aggregate["backbone"].eq("kotha2020_rjb")
        & aggregate["minimum_mw"].eq(4.5)
    ].sort_values("period_s")
    supported = terms[
        terms["backbone"].eq("bindi2014_rhyp")
        & terms["minimum_mw"].eq(4.5)
        & terms["start_date"].eq("2015-01-01")
        & terms["n_records"].ge(15)
    ]
    transfer3 = transfer_metrics[transfer_metrics["period_s"].eq(3.0)].iloc[0]
    spatial3 = spatial_metrics[
        spatial_metrics["period_s"].eq(3.0)
        & spatial_metrics["scope"].isin(["overall", "event_fold_mean"])
    ].sort_values(["validation_design", "model"])
    lines = [
        "# ESM external station-field validation",
        "",
        "## Data boundary",
        "",
        "- Source: Engineering Strong-Motion Database (ESM) v3 flatfile web service.",
        f"- Local snapshot SHA-256: `{metadata['sha256']}`.",
        f"- Snapshot rows: {len(raw):,}; events through {metadata['maximum_event_time']}.",
        "- The primary sample retains European-Mediterranean shallow events with 4.5 <= Mw <= 7.6 from 2015 onward, after publication of the fixed Bindi et al. (2014) backbone.",
        "- As-recorded horizontal geometric-mean spectra match the primary backbone component. Kotha et al. (2020), its Mw <= 7.4 range and RotD50 are retained as a smaller response-component and backbone sensitivity.",
        "- Event-station graph components are solved separately. Independent offsets are removed within intersecting train-test component pairs before pooled comparison.",
        "",
        "## Primary post-2015 replication",
        "",
        f"- Supported field: {supported['station_key'].nunique():,} stations with at least 15 records across {supported['period_s'].nunique()} periods.",
    ]
    for row in primary.itertuples():
        lines.append(
            f"- SA({row.period_s:g} s): mean held-event correlation {row.pearson:.3f} "
            f"(station-bootstrap 95% CI {row.pearson_ci_low:.3f} to {row.pearson_ci_high:.3f}); "
            f"RMSE gain {row.rmse_gain_pct:.1f}% "
            f"({row.rmse_gain_pct_ci_low:.1f}% to {row.rmse_gain_pct_ci_high:.1f}%)."
        )
    lines.extend(["", "## Sensitivities", ""])
    for row in sensitivity.itertuples():
        lines.append(
            f"- Kotha/RotD50 SA({row.period_s:g} s): correlation {row.pearson:.3f}; "
            f"RMSE gain {row.rmse_gain_pct:.1f}%."
        )
    lines.extend(
        [
            "",
            "## Cross-region transfer boundary",
            "",
            f"- A frozen Japan K-NET model using only log10(Vs30) and elevation gives SA(3.0 s) "
            f"correlation {transfer3['pearson']:.3f} and RMSE gain {transfer3['rmse_gain_pct']:.1f}% in Europe.",
            "- The adverse transfer is retained. Japan-specific basin-depth and volcanic-front variables are unavailable in ESM, so this test evaluates only the two harmonised variables and does not invalidate the Japan regional model.",
            "- The external result supports repeatability of long-period station fields across regions. It does not support a globally transferable station predictor.",
            "",
            "## ESM local spatial prediction",
            "",
            f"- K-means requested {REQUESTED_SPATIAL_BLOCKS} geographic blocks. The fixed minimum of "
            f"{MINIMUM_SPATIAL_BLOCK_STATIONS} stations per block selected "
            f"{int(spatial3['n_spatial_blocks'].iloc[0])} blocks without using response values.",
            "- Three feature sets were fixed before evaluation: harmonised VS30 and elevation; all available ESM site metadata; and the same metadata plus coordinates. Hyperparameters match the Japan analysis and are not tuned on ESM outcomes.",
            "- The strict design estimates the training field from four event groups, predicts spatially held stations, and scores predictions against the fifth event group. Training targets are centred only on spatial-training stations; out-of-fold predictions are centred without response labels.",
        ]
    )
    for row in spatial3.itertuples():
        lines.append(
            f"- {row.validation_design}, {row.model}, SA(3.0 s): correlation {row.pearson:.3f} "
            f"(95% CI {row.pearson_ci_low:.3f} to {row.pearson_ci_high:.3f}); RMSE gain "
            f"{row.rmse_gain_pct:.1f}% ({row.rmse_gain_pct_ci_low:.1f}% to "
            f"{row.rmse_gain_pct_ci_high:.1f}%)."
        )
    lines.extend(
        [
            "- All spatial results, including adverse folds and intervals, are retained. They test geographic extrapolation of available public ESM variables and do not test unavailable basin-depth variables.",
            "",
            "## Selection table",
            "",
            "```text",
            selection.to_csv(index=False, sep="\t").rstrip(),
            "```",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_outputs(selection: pd.DataFrame, terms: pd.DataFrame,
                     repeatability_metrics: pd.DataFrame, transfer_metrics: pd.DataFrame,
                     spatial_metrics: pd.DataFrame) -> None:
    primary = repeatability_metrics[
        repeatability_metrics["scope"].eq("fold_mean")
        & repeatability_metrics["backbone"].eq("bindi2014_rhyp")
        & repeatability_metrics["minimum_mw"].eq(4.5)
        & repeatability_metrics["start_date"].eq("2015-01-01")
    ]
    assert set(primary["period_s"]) == set(PERIODS)
    assert primary["pearson"].notna().all() and primary["rmse_gain_pct"].notna().all()
    assert len(terms) > 0 and terms["station_key"].notna().all()
    assert selection["n_records"].ge(0).all()
    assert set(transfer_metrics["period_s"]) == set(PERIODS)
    spatial_overall = spatial_metrics[spatial_metrics["scope"].eq("overall")]
    strict_overall = spatial_metrics[spatial_metrics["scope"].eq("event_fold_mean")]
    assert set(spatial_overall["period_s"]) == set(PERIODS)
    assert set(strict_overall["period_s"]) == set(PERIODS)
    assert set(spatial_overall["model"]) == set(MODEL_FEATURES)
    assert set(strict_overall["model"]) == set(MODEL_FEATURES)
    assert spatial_metrics["n_spatial_blocks"].ge(2).all()
    for path in [OUT_SELECTION, OUT_DECOMPOSITION, OUT_TERMS, OUT_REPEATABILITY, OUT_REPEATABILITY_PAIRS,
                 OUT_FIELD_COMPARISON, OUT_TRANSFER_METRICS, OUT_TRANSFER_PREDICTIONS,
                 OUT_SPATIAL_METRICS, OUT_SPATIAL_PREDICTIONS,
                 OUT_FIGURE_PDF, OUT_FIGURE_PNG, OUT_AUDIT]:
        assert path.exists() and path.stat().st_size > 0, path


def run(args: argparse.Namespace) -> None:
    download_data(args.data, args.force_download)
    raw = pd.read_csv(args.data, sep=";", low_memory=False)
    missing = sorted(set(QUERY_FIELDS) - set(raw.columns))
    if missing:
        raise ValueError(f"ESM query is missing fields: {missing}")
    metadata = {
        "source_url": "https://esm-db.eu/esmws/flatfile/1/",
        "query_url": QUERY_URL,
        "sha256": sha256(args.data),
        "n_rows": len(raw),
        "maximum_event_time": str(pd.to_datetime(raw["event_time"], errors="coerce").max()),
    }
    DATA_META.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    analyses = [
        ("bindi2014_rhyp", 4.5, "2015-01-01", 15, 5, 10),
        ("bindi2014_rhyp", 5.0, "2015-01-01", 15, 5, 10),
        ("kotha2020_rjb", 4.5, None, 10, 3, 3),
    ]
    selection_rows = []
    term_blocks = []
    decomposition_metric_blocks = []
    repeatability_metric_blocks = []
    repeatability_pair_blocks = []
    for backbone, minimum_mw, start_date, min_train, min_test, min_pair in analyses:
        selected, audit = select_records(raw, backbone, minimum_mw, start_date)
        selection_rows.extend(audit)
        frames = residual_frames(selected, backbone)
        terms, decomposition_metrics = full_station_terms(frames, backbone, minimum_mw, start_date)
        metrics, pairs = repeatability(
            frames, backbone, minimum_mw, start_date, min_train, min_test, min_pair, args.seed
        )
        term_blocks.append(terms)
        decomposition_metric_blocks.append(decomposition_metrics)
        repeatability_metric_blocks.append(metrics)
        repeatability_pair_blocks.append(pairs)
        print(
            f"ESM analysis: {backbone} Mw={minimum_mw:g}-{MAXIMUM_MW[backbone]:g} "
            f"start={start_date or 'all'} "
            f"records={len(selected):,}", flush=True,
        )

    selection = pd.DataFrame(selection_rows)
    terms = pd.concat(term_blocks, ignore_index=True)
    decomposition_metrics = pd.concat(decomposition_metric_blocks, ignore_index=True)
    fold_metrics = pd.concat(repeatability_metric_blocks, ignore_index=True)
    pairs = pd.concat(repeatability_pair_blocks, ignore_index=True)
    repeatability_metrics = aggregate_repeatability(
        fold_metrics, pairs, args.bootstrap_replicates, args.seed
    )
    comparison, japan = field_comparison(terms, minimum_records=15)
    transfer_metrics, transfer_predictions = cross_region_transfer(
        terms, japan, args.bootstrap_replicates, args.seed
    )
    spatial_metrics, spatial_predictions = local_spatial_prediction(
        terms, pairs, args.bootstrap_replicates, args.seed
    )

    selection.to_csv(OUT_SELECTION, index=False)
    decomposition_metrics.to_csv(OUT_DECOMPOSITION, index=False)
    terms.to_csv(OUT_TERMS, index=False)
    repeatability_metrics.to_csv(OUT_REPEATABILITY, index=False)
    pairs.to_csv(OUT_REPEATABILITY_PAIRS, index=False)
    comparison.to_csv(OUT_FIELD_COMPARISON, index=False)
    transfer_metrics.to_csv(OUT_TRANSFER_METRICS, index=False)
    transfer_predictions.to_csv(OUT_TRANSFER_PREDICTIONS, index=False)
    spatial_metrics.to_csv(OUT_SPATIAL_METRICS, index=False)
    spatial_predictions.to_csv(OUT_SPATIAL_PREDICTIONS, index=False)
    plot_results(terms, repeatability_metrics, pairs, comparison)
    write_audit(
        raw, metadata, selection, terms, repeatability_metrics, transfer_metrics, spatial_metrics
    )
    validate_outputs(selection, terms, repeatability_metrics, transfer_metrics, spatial_metrics)

    primary = repeatability_metrics[
        repeatability_metrics["scope"].eq("fold_mean")
        & repeatability_metrics["backbone"].eq("bindi2014_rhyp")
        & repeatability_metrics["minimum_mw"].eq(4.5)
        & repeatability_metrics["start_date"].eq("2015-01-01")
    ].sort_values("period_s")
    print(primary[["period_s", "pearson", "rmse_gain_pct"]].to_string(index=False), flush=True)
    print(f"wrote {OUT_AUDIT}", flush=True)
    print(f"wrote {OUT_FIGURE_PDF}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DATA_FILE)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
