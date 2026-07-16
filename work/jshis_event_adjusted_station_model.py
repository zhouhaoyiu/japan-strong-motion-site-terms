#!/usr/bin/env python3
"""Event-adjusted, direct-period station residual analysis for the CEE paper.

The analysis uses the public J-SHIS/NIED strong-motion flatfile and the MF2013
coefficient table.  For each official response-spectrum period it:

1. computes RotD100 residuals against the MF2013 backbone and D1400/AVS30 site terms;
2. separates the global intercept, event term, and station term with a two-way
   additive fixed-effect decomposition;
3. predicts the zero-centred station term in spatially blocked cross-validation;
4. propagates the out-of-fold prediction to matched official response-spectrum
   ordinates after converting the Vs=400 m/s engineering-bedrock ordinate to
   the station AVS30 reference condition.

The AI term is evaluated with the public rule-based geographic approximation as
a sensitivity analysis.  PH is event-constant at a given period and is absorbed
by the event term; it is not used to interpret the global MF2013 intercept.
"""

from __future__ import annotations

import argparse
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(__file__).resolve().parents[1]
ARTICLE_DIR = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT_DIR = ARTICLE_DIR / "supplement"
FIGURE_DIR = ARTICLE_DIR / "figures"

DEFAULT_FLATFILE = ROOT / "work" / "external_data" / "jshis_gmf" / "flatfile_sub1-v2024.zip"
DEFAULT_COEFFICIENTS = ROOT / "work" / "external_data" / "jshis_mf2013" / "MF13rev_coefs.csv"
DEFAULT_RESPONSE_MAP = ROOT / "work" / "external_data" / "jshis_respmap" / "P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip"

OUT_STATION_TERMS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_terms.csv"
OUT_DECOMPOSITION = SUPPLEMENT_DIR / "jshis_event_station_decomposition_metrics.csv"
OUT_MODEL_PREDICTIONS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv"
OUT_MODEL_METRICS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_metrics.csv"
OUT_EVENT_HOLDOUT = SUPPLEMENT_DIR / "jshis_event_holdout_station_repeatability.csv"
OUT_HAZARD_VALUES = SUPPLEMENT_DIR / "jshis_event_adjusted_surface_spectrum_values.csv"
OUT_HAZARD_SUMMARY = SUPPLEMENT_DIR / "jshis_event_adjusted_surface_spectrum_summary.csv"
OUT_CITY_CASES = SUPPLEMENT_DIR / "jshis_event_adjusted_city_nearest_cases.csv"
OUT_AUDIT = SUPPLEMENT_DIR / "jshis_event_adjusted_station_analysis.md"
OUT_COMPONENT_SUMMARY = SUPPLEMENT_DIR / "jshis_response_component_audit.csv"
OUT_COMPONENT_STATIONS = SUPPLEMENT_DIR / "jshis_response_component_station_terms.csv"
OUT_COMPONENT_AUDIT = SUPPLEMENT_DIR / "jshis_response_component_audit.md"
OUT_FIGURE_PDF = FIGURE_DIR / "figure_event_adjusted_multiperiod.pdf"
OUT_FIGURE_PNG = FIGURE_DIR / "figure_event_adjusted_multiperiod.png"
OUT_OVERVIEW_PDF = FIGURE_DIR / "figure_event_adjusted_overview.pdf"
OUT_OVERVIEW_PNG = FIGURE_DIR / "figure_event_adjusted_overview.png"
OUT_STRUCTURE_PDF = FIGURE_DIR / "figure_event_adjusted_site_structure.pdf"
OUT_STRUCTURE_PNG = FIGURE_DIR / "figure_event_adjusted_site_structure.png"
OUT_CITY_PDF = FIGURE_DIR / "figure_event_adjusted_city_cases.pdf"
OUT_CITY_PNG = FIGURE_DIR / "figure_event_adjusted_city_cases.png"

G_IN_CM_S2 = 980.665


@dataclass(frozen=True)
class PeriodSpec:
    code: str
    period_s: float
    coeff_key: str
    rotd50_col: str
    rotd100_col: str

    @property
    def label(self) -> str:
        return f"SA({self.period_s:g} s)"


PERIODS = [
    PeriodSpec("P0010", 0.1, "0.10", "rsaccrd050d005t0010", "rsaccrd100d005t0010"),
    PeriodSpec("P0020", 0.2, "0.20", "rsaccrd050d005t0020", "rsaccrd100d005t0020"),
    PeriodSpec("P0030", 0.3, "0.30", "rsaccrd050d005t0030", "rsaccrd100d005t0030"),
    PeriodSpec("P0050", 0.5, "0.50", "rsaccrd050d005t0050", "rsaccrd100d005t0050"),
    PeriodSpec("P0100", 1.0, "1.00", "rsaccrd050d005t0100", "rsaccrd100d005t0100"),
    PeriodSpec("P0200", 2.0, "2.00", "rsaccrd050d005t0200", "rsaccrd100d005t0200"),
    PeriodSpec("P0300", 3.0, "3.00", "rsaccrd050d005t0300", "rsaccrd100d005t0300"),
    PeriodSpec("P0500", 5.0, "5.00", "rsaccrd050d005t0500", "rsaccrd100d005t0500"),
]

PROBABILITY_COLUMNS = {
    "T50_P02_SA": "50y_2pct",
    "T50_P05_SA": "50y_5pct",
    "T50_P10_SA": "50y_10pct",
    "T50_P39_SA": "50y_39pct",
}

SITE_COLUMNS = [
    "siteid2",
    "site_code",
    "lon",
    "lat",
    "elevation",
    "sensor_depth_glminus",
    "obs_network_id",
    "installation_situation_id",
    "dist_vf_mf13_nejapan",
    "dist_vf_mf13_swjapan",
    "vs10",
    "vs20",
    "vs30",
    "avs30",
    "meshcode250",
    "meshcode3",
    "d1100",
    "d1400",
    "d1700",
    "d2100",
    "dbase",
]

SOURCE_COLUMNS = [
    "eq_source_id",
    "segment_idx",
    "jem_lat",
    "jem_lon",
    "jem_depth",
    "mw",
    "rake1",
    "eq_location_type_id",
]

PHYSICAL_FEATURES = [
    "log_vs10",
    "log_vs20",
    "log_vs30",
    "log_avs30",
    "log_d1100",
    "log_d1400",
    "log_d1700",
    "log_d2100",
    "log_dbase",
    "elevation",
    "sensor_depth_glminus",
]

SPATIAL_FEATURES = [
    "lon",
    "lat",
    "dist_vf_mf13_nejapan",
    "dist_vf_mf13_swjapan",
]

CATEGORICAL_FEATURES = ["network_label"]

PRIMARY_INSTALLATION_SITUATION_ID = 1  # GROUND in the official flatfile code table.

CITIES = {
    "Tokyo": (139.6917, 35.6895),
    "Osaka": (135.5023, 34.6937),
    "Sapporo": (141.3545, 43.0618),
    "Sendai": (140.8719, 38.2682),
    "Fukuoka": (130.4017, 33.5904),
}


def normalize_mesh(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def load_coefficients(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame.columns = [str(col).strip() for col in frame.columns]
    frame["Period"] = frame["Period"].astype(str).str.strip()
    for col in frame.columns:
        if col != "Period":
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame.set_index("Period")


def load_metadata(flatfile: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with zipfile.ZipFile(flatfile) as archive:
        site = pd.read_csv(archive.open("site_schema.tsv"), sep="\t", usecols=SITE_COLUMNS)
        source = pd.read_csv(archive.open("source_schema.tsv"), sep="\t", usecols=SOURCE_COLUMNS)
    for col in SITE_COLUMNS:
        if col != "site_code":
            site[col] = pd.to_numeric(site[col], errors="coerce")
    site["siteid2"] = site["siteid2"].astype("Int64")
    site = site.sort_values(["siteid2", "site_code"]).drop_duplicates("siteid2", keep="last")
    site["network_label"] = site["obs_network_id"].map({1: "K-NET", 2: "KiK-net"}).fillna("other")
    site["meshcode250_str"] = site["meshcode250"].map(normalize_mesh)
    for col in ["vs10", "vs20", "vs30", "avs30", "d1100", "d1400", "d1700", "d2100", "dbase"]:
        site[f"log_{col}"] = np.log10(pd.to_numeric(site[col], errors="coerce").where(lambda x: x > 0))

    for col in SOURCE_COLUMNS:
        source[col] = pd.to_numeric(source[col], errors="coerce")
    source["eq_source_id"] = source["eq_source_id"].astype("Int64")
    source = source.sort_values(["eq_source_id", "segment_idx"]).drop_duplicates("eq_source_id", keep="first")
    return site, source


def source_class_terms(source_class: pd.Series, coeff: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    kinds = pd.to_numeric(source_class, errors="coerce")
    b = np.select(
        [kinds.eq(1), kinds.eq(2), kinds.eq(3)],
        [coeff["b1"], coeff["b2"], coeff["b3"]],
        default=np.nan,
    )
    c = np.select(
        [kinds.eq(1), kinds.eq(2), kinds.eq(3)],
        [coeff["c1"], coeff["c2"], coeff["c3"]],
        default=np.nan,
    )
    return b, c


def mf2013_predictions(frame: pd.DataFrame, coeff: pd.Series) -> tuple[pd.Series, pd.Series, pd.Series]:
    mw = np.minimum(pd.to_numeric(frame["mw"], errors="coerce"), 8.2)
    distance = np.maximum(pd.to_numeric(frame["fault_dist"], errors="coerce"), 1.0)
    b, c = source_class_terms(frame["eq_location_type_id"], coeff)
    basic = (
        coeff["a"] * (mw - 16.0) ** 2
        + b * distance
        + c
        - np.log10(distance + coeff["d"] * 10.0 ** (0.5 * mw))
    )
    basic = pd.Series(basic, index=frame.index, dtype=float)

    d1400 = pd.to_numeric(frame["d1400"], errors="coerce")
    avs30 = pd.to_numeric(frame["avs30"], errors="coerce")
    gd = coeff["pd"] * np.log10(np.maximum(coeff["Dlmin"], d1400) / 300.0)
    gs = coeff["ps"] * np.log10(np.minimum(coeff["Vsmax"], avs30.where(avs30 > 0)) / 350.0)
    site = basic + gd + gs

    depth = pd.to_numeric(frame["jem_depth"], errors="coerce")
    source_lon = pd.to_numeric(frame["jem_lon"], errors="coerce")
    site_lon = pd.to_numeric(frame["lon"], errors="coerce")
    site_lat = pd.to_numeric(frame["lat"], errors="coerce")
    kind = pd.to_numeric(frame["eq_location_type_id"], errors="coerce")
    subduction = kind.isin([2, 3])
    ne = subduction & depth.gt(30.0) & site_lat.ge(36.0) & source_lon.ge(136.9)
    sw = subduction & depth.ge(60.0) & site_lon.lt(136.9) & source_lon.lt(136.9)
    ai = pd.Series(0.0, index=frame.index, dtype=float)
    ne_distance = pd.to_numeric(frame["dist_vf_mf13_nejapan"], errors="coerce")
    sw_distance = pd.to_numeric(frame["dist_vf_mf13_swjapan"], errors="coerce").clip(upper=75.0)
    ai.loc[ne] = coeff["gNE"] * ne_distance.loc[ne] * (depth.loc[ne] - 30.0)
    ai.loc[sw] = coeff["gSW"] * sw_distance.loc[sw] * (depth.loc[sw] - 30.0)
    return basic, site, site + ai


def load_record_residuals(
    flatfile: Path,
    coefficients: pd.DataFrame,
    site: pd.DataFrame,
    source: pd.DataFrame,
    chunksize: int,
    installation_situation_id: int | None = PRIMARY_INSTALLATION_SITUATION_ID,
) -> dict[float, pd.DataFrame]:
    response_cols = [column for spec in PERIODS for column in (spec.rotd50_col, spec.rotd100_col)]
    usecols = ["smrec_id", "site_id", "eq_source_id", "fault_dist", *response_cols]
    buffers: dict[float, list[pd.DataFrame]] = {spec.period_s: [] for spec in PERIODS}
    joined_rows = 0
    with zipfile.ZipFile(flatfile) as archive:
        reader = pd.read_csv(
            archive.open("smrec_schema.tsv"),
            sep="\t",
            usecols=usecols,
            chunksize=chunksize,
            low_memory=False,
        )
        for chunk_number, chunk in enumerate(reader, start=1):
            chunk["site_id"] = pd.to_numeric(chunk["site_id"], errors="coerce").astype("Int64")
            chunk["siteid2"] = (chunk["site_id"] // 10).astype("Int64")
            chunk["eq_source_id"] = pd.to_numeric(chunk["eq_source_id"], errors="coerce").astype("Int64")
            joined = chunk.merge(site, on="siteid2", how="left", validate="many_to_one")
            joined = joined.merge(source, on="eq_source_id", how="left", validate="many_to_one")
            joined_rows += len(joined)
            for spec in PERIODS:
                observed_raw = pd.to_numeric(joined[spec.rotd100_col], errors="coerce")
                observed = np.log10(observed_raw.where(observed_raw > 0))
                observed_rotd50_raw = pd.to_numeric(joined[spec.rotd50_col], errors="coerce")
                observed_rotd50 = np.log10(observed_rotd50_raw.where(observed_rotd50_raw > 0))
                basic, site_pred, site_ai = mf2013_predictions(joined, coefficients.loc[spec.coeff_key])
                usable = (
                    observed.notna()
                    & basic.notna()
                    & site_pred.notna()
                    & site_ai.notna()
                    & np.isfinite(observed)
                    & np.isfinite(basic)
                    & np.isfinite(site_pred)
                    & np.isfinite(site_ai)
                    & pd.to_numeric(joined["fault_dist"], errors="coerce").gt(0)
                    & pd.to_numeric(joined["avs30"], errors="coerce").gt(0)
                    & joined["eq_location_type_id"].isin([1, 2, 3])
                )
                if installation_situation_id is not None:
                    usable &= joined["installation_situation_id"].eq(installation_situation_id)
                block = joined.loc[usable, ["smrec_id", "siteid2", "eq_source_id", "fault_dist"]].copy()
                block["residual_basic"] = observed.loc[usable].to_numpy() - basic.loc[usable].to_numpy()
                block["residual_site"] = observed.loc[usable].to_numpy() - site_pred.loc[usable].to_numpy()
                block["residual_site_ai"] = observed.loc[usable].to_numpy() - site_ai.loc[usable].to_numpy()
                block["residual_site_rotd50"] = observed_rotd50.loc[usable].to_numpy() - site_pred.loc[usable].to_numpy()
                buffers[spec.period_s].append(block)
            print(f"records: chunk={chunk_number} joined={joined_rows:,}", flush=True)
    return {period: pd.concat(parts, ignore_index=True) for period, parts in buffers.items()}


def weighted_rmse(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.average((y_true - y_pred) ** 2, weights=weights)))


def weighted_mae(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(np.abs(y_true - y_pred), weights=weights))


def two_way_decomposition(
    frame: pd.DataFrame,
    residual_col: str,
    tolerance: float = 1.0e-10,
    max_iterations: int = 5_000,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    clean = frame[["siteid2", "eq_source_id", residual_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    station_codes, station_ids = pd.factorize(clean["siteid2"], sort=True)
    event_codes, event_ids = pd.factorize(clean["eq_source_id"], sort=True)
    values = clean[residual_col].to_numpy(float)
    station_counts = np.bincount(station_codes).astype(float)
    event_counts = np.bincount(event_codes).astype(float)
    station_effect = np.zeros(len(station_counts), dtype=float)
    event_effect = np.zeros(len(event_counts), dtype=float)
    intercept = float(np.mean(values))

    max_change = math.inf
    for iteration in range(1, max_iterations + 1):
        old_station = station_effect.copy()
        old_event = event_effect.copy()
        station_effect = np.bincount(
            station_codes,
            weights=values - intercept - event_effect[event_codes],
            minlength=len(station_counts),
        ) / station_counts
        shift = float(np.average(station_effect, weights=station_counts))
        station_effect -= shift
        intercept += shift

        event_effect = np.bincount(
            event_codes,
            weights=values - intercept - station_effect[station_codes],
            minlength=len(event_counts),
        ) / event_counts
        shift = float(np.average(event_effect, weights=event_counts))
        event_effect -= shift
        intercept += shift
        max_change = max(
            float(np.max(np.abs(station_effect - old_station))),
            float(np.max(np.abs(event_effect - old_event))),
        )
        if max_change < tolerance:
            break

    fitted = intercept + station_effect[station_codes] + event_effect[event_codes]
    remainder = values - fitted
    station_table = pd.DataFrame(
        {
            "siteid2": pd.Series(station_ids, dtype="Int64"),
            "n_records": station_counts.astype(int),
            "station_effect_log10": station_effect,
        }
    )
    event_table = pd.DataFrame(
        {
            "eq_source_id": pd.Series(event_ids, dtype="Int64"),
            "n_records": event_counts.astype(int),
            "event_effect_log10": event_effect,
        }
    )
    metrics = {
        "n_records": float(len(clean)),
        "n_stations": float(len(station_table)),
        "n_events": float(len(event_table)),
        "global_intercept_log10": intercept,
        "raw_mae_log10": float(mean_absolute_error(np.zeros_like(values), values)),
        "raw_rmse_log10": float(np.sqrt(mean_squared_error(np.zeros_like(values), values))),
        "remainder_mae_log10": float(mean_absolute_error(np.zeros_like(remainder), remainder)),
        "remainder_rmse_log10": float(np.sqrt(mean_squared_error(np.zeros_like(remainder), remainder))),
        "station_effect_weighted_std_log10": float(np.sqrt(np.average(station_effect**2, weights=station_counts))),
        "event_effect_weighted_std_log10": float(np.sqrt(np.average(event_effect**2, weights=event_counts))),
        "iterations": float(iteration),
        "max_parameter_change": max_change,
        "weighted_station_mean": float(np.average(station_effect, weights=station_counts)),
        "weighted_event_mean": float(np.average(event_effect, weights=event_counts)),
    }
    return station_table, event_table, metrics


def decompose_all_periods(
    records: dict[float, pd.DataFrame],
    site: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    term_rows = []
    metric_rows = []
    model_columns = {
        "mf2013_basic": "residual_basic",
        "mf2013_site": "residual_site",
        "mf2013_site_ai_sensitivity": "residual_site_ai",
    }
    site_metadata = site.drop_duplicates("siteid2")
    for spec in PERIODS:
        frame = records[spec.period_s]
        for model_name, residual_col in model_columns.items():
            station_terms, _, metrics = two_way_decomposition(frame, residual_col)
            station_terms["period_s"] = spec.period_s
            station_terms["period_code"] = spec.code
            station_terms["model"] = model_name
            station_terms["response_component"] = "RotD100"
            term_rows.append(station_terms)
            metric_rows.append(
                {
                    "period_s": spec.period_s,
                    "period_code": spec.code,
                    "model": model_name,
                    "response_component": "RotD100",
                    **metrics,
                }
            )
        print(f"decomposition: period={spec.period_s:g}s", flush=True)
    terms = pd.concat(term_rows, ignore_index=True).merge(site_metadata, on="siteid2", how="left", validate="many_to_one")
    metrics = pd.DataFrame(metric_rows)

    correlation_rows = []
    supported = terms[terms["n_records"].ge(20)].copy()
    for (period_s, model), sub in supported.groupby(["period_s", "model"]):
        for feature in ["avs30", "d1400", "dbase"]:
            valid = pd.to_numeric(sub[feature], errors="coerce").gt(0) & sub["station_effect_log10"].notna()
            rho = sub.loc[valid, [feature, "station_effect_log10"]].corr(method="spearman").iloc[0, 1]
            correlation_rows.append(
                {
                    "period_s": period_s,
                    "model": model,
                    "association_feature": feature,
                    "association_n_stations": int(valid.sum()),
                    "station_effect_spearman": float(rho),
                }
            )
    correlations = pd.DataFrame(correlation_rows)
    metrics = metrics.merge(correlations, on=["period_s", "model"], how="left")
    return terms, metrics


def audit_response_components(
    records: dict[float, pd.DataFrame],
    station_terms: pd.DataFrame,
    decomposition: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    station_rows = []
    for spec in PERIODS:
        frame = records[spec.period_s]
        rotd50_terms, _, rotd50_metrics = two_way_decomposition(frame, "residual_site_rotd50")
        rotd100_terms = station_terms[
            station_terms["period_s"].eq(spec.period_s) & station_terms["model"].eq("mf2013_site")
        ][["siteid2", "n_records", "station_effect_log10"]].rename(
            columns={
                "n_records": "n_records_rotd100",
                "station_effect_log10": "station_effect_rotd100_log10",
            }
        )
        rotd50_terms = rotd50_terms.rename(
            columns={
                "n_records": "n_records_rotd50",
                "station_effect_log10": "station_effect_rotd50_log10",
            }
        )
        paired = rotd100_terms.merge(rotd50_terms, on="siteid2", how="inner", validate="one_to_one")
        paired["period_s"] = spec.period_s
        paired["period_code"] = spec.code
        paired["station_difference_rotd50_minus_rotd100_log10"] = (
            paired["station_effect_rotd50_log10"] - paired["station_effect_rotd100_log10"]
        )
        paired["weight_records"] = paired[["n_records_rotd50", "n_records_rotd100"]].min(axis=1)
        station_rows.append(paired)

        differences = paired["station_difference_rotd50_minus_rotd100_log10"].to_numpy(float)
        weights = paired["weight_records"].to_numpy(float)
        record_ratio = np.power(
            10.0,
            frame["residual_site_rotd50"].to_numpy(float) - frame["residual_site"].to_numpy(float),
        )
        rotd100_metrics = decomposition[
            decomposition["period_s"].eq(spec.period_s) & decomposition["model"].eq("mf2013_site")
        ].iloc[0]
        summary_rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "n_records_rotd100": int(rotd100_metrics["n_records"]),
                "n_records_rotd50": int(rotd50_metrics["n_records"]),
                "n_paired_stations": len(paired),
                "station_pearson": paired[
                    ["station_effect_rotd100_log10", "station_effect_rotd50_log10"]
                ].corr(method="pearson").iloc[0, 1],
                "station_spearman": paired[
                    ["station_effect_rotd100_log10", "station_effect_rotd50_log10"]
                ].corr(method="spearman").iloc[0, 1],
                "station_difference_weighted_mean_log10": float(np.average(differences, weights=weights)),
                "station_difference_weighted_rmse_log10": weighted_rmse(
                    paired["station_effect_rotd100_log10"].to_numpy(float),
                    paired["station_effect_rotd50_log10"].to_numpy(float),
                    weights,
                ),
                "station_difference_abs_q95_log10": float(np.quantile(np.abs(differences), 0.95)),
                "record_rotd50_over_rotd100_q05": float(np.nanquantile(record_ratio, 0.05)),
                "record_rotd50_over_rotd100_q50": float(np.nanquantile(record_ratio, 0.50)),
                "record_rotd50_over_rotd100_q95": float(np.nanquantile(record_ratio, 0.95)),
                "remainder_rmse_rotd100_log10": float(rotd100_metrics["remainder_rmse_log10"]),
                "remainder_rmse_rotd50_log10": float(rotd50_metrics["remainder_rmse_log10"]),
            }
        )
        print(f"response-component audit: period={spec.period_s:g}s", flush=True)

    summary = pd.DataFrame(summary_rows)
    stations = pd.concat(station_rows, ignore_index=True)
    summary.to_csv(OUT_COMPONENT_SUMMARY, index=False)
    stations.to_csv(OUT_COMPONENT_STATIONS, index=False)
    sa3 = summary[summary["period_s"].eq(3.0)].iloc[0]
    lines = [
        "# Response-component audit",
        "",
        "MF2013 predicts the maximum norm of the two horizontal oscillator responses. This coordinate is equivalent to RotD100, so RotD100 is used in the primary residual analysis. RotD50 is retained as a response-component sensitivity test and as the observed coordinate for the geometric-mean Zhao model comparison.",
        "",
        "## SA(3.0 s)",
        "",
        f"- Paired station-term Pearson correlation: {sa3['station_pearson']:.4f}.",
        f"- Record-weighted station-term difference RMSE: {sa3['station_difference_weighted_rmse_log10']:.4f} log10 units.",
        f"- Station-term absolute-difference 95th percentile: {sa3['station_difference_abs_q95_log10']:.4f} log10 units.",
        f"- Record-level RotD50/RotD100 median: {sa3['record_rotd50_over_rotd100_q50']:.4f}.",
        f"- RotD100 and RotD50 remainder RMSE: {sa3['remainder_rmse_rotd100_log10']:.4f} and {sa3['remainder_rmse_rotd50_log10']:.4f} log10 units.",
        "",
        "The response-component choice changes the absolute spectral coordinate and residual dispersion, while the zero-centred station field is stable.",
    ]
    OUT_COMPONENT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary, stations


def event_holdout_repeatability(
    records: dict[float, pd.DataFrame],
    n_splits: int,
    seed: int,
    min_train_records: int = 15,
    min_test_records: int = 5,
) -> pd.DataFrame:
    """Measure whether two-way station effects recur in disjoint event sets."""

    rows = []
    for spec in PERIODS:
        frame = records[spec.period_s][["siteid2", "eq_source_id", "residual_site"]].copy()
        events = np.sort(frame["eq_source_id"].dropna().unique())
        # The same event partition is used at every period and for each backbone.
        rng = np.random.default_rng(seed)
        shuffled = rng.permutation(events)
        event_folds = {event_id: index % n_splits for index, event_id in enumerate(shuffled)}
        fold_ids = frame["eq_source_id"].map(event_folds).to_numpy(int)
        period_rows = []
        for fold in range(n_splits):
            train = frame.loc[fold_ids != fold].copy()
            test = frame.loc[fold_ids == fold].copy()
            train_terms, _, train_fit = two_way_decomposition(train, "residual_site")
            test_terms, _, test_fit = two_way_decomposition(test, "residual_site")
            train_terms = train_terms.rename(
                columns={
                    "station_effect_log10": "train_station_term",
                    "n_records": "n_train_records",
                }
            )
            test_terms = test_terms.rename(
                columns={
                    "station_effect_log10": "test_station_term",
                    "n_records": "n_test_records",
                }
            )
            paired = train_terms.merge(test_terms, on="siteid2", how="inner", validate="one_to_one")
            paired = paired[
                paired["n_train_records"].ge(min_train_records)
                & paired["n_test_records"].ge(min_test_records)
            ].copy()
            train_values = paired["train_station_term"].to_numpy(float, copy=True)
            test_values = paired["test_station_term"].to_numpy(float, copy=True)
            weights = paired["n_test_records"].to_numpy(float)
            train_values -= np.average(train_values, weights=paired["n_train_records"].to_numpy(float))
            test_values -= np.average(test_values, weights=weights)
            baseline_rmse = weighted_rmse(test_values, np.zeros_like(test_values), weights)
            prediction_rmse = weighted_rmse(test_values, train_values, weights)
            row = {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "scope": "fold",
                "fold": fold,
                "decomposition_method": "two_way_fixed_effects",
                "event_split_seed": seed,
                "n_train_events": int(train["eq_source_id"].nunique()),
                "n_test_events": int(test["eq_source_id"].nunique()),
                "n_paired_stations": len(paired),
                "test_record_weight": int(weights.sum()),
                "train_iterations": int(train_fit["iterations"]),
                "test_iterations": int(test_fit["iterations"]),
                "train_max_parameter_change": float(train_fit["max_parameter_change"]),
                "test_max_parameter_change": float(test_fit["max_parameter_change"]),
                "train_test_station_correlation": float(np.corrcoef(train_values, test_values)[0, 1]),
                "zero_baseline_rmse_log10": baseline_rmse,
                "train_term_prediction_rmse_log10": prediction_rmse,
                "rmse_reduction_vs_zero_pct": 100.0 * (baseline_rmse - prediction_rmse) / baseline_rmse,
            }
            rows.append(row)
            period_rows.append(row)
        period_frame = pd.DataFrame(period_rows)
        rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "scope": "fold_mean",
                "fold": -1,
                "decomposition_method": "two_way_fixed_effects",
                "event_split_seed": seed,
                "n_train_events": float(period_frame["n_train_events"].mean()),
                "n_test_events": float(period_frame["n_test_events"].mean()),
                "n_paired_stations": float(period_frame["n_paired_stations"].mean()),
                "test_record_weight": float(period_frame["test_record_weight"].sum()),
                "train_iterations": int(period_frame["train_iterations"].max()),
                "test_iterations": int(period_frame["test_iterations"].max()),
                "train_max_parameter_change": float(period_frame["train_max_parameter_change"].max()),
                "test_max_parameter_change": float(period_frame["test_max_parameter_change"].max()),
                "train_test_station_correlation": float(period_frame["train_test_station_correlation"].mean()),
                "zero_baseline_rmse_log10": float(period_frame["zero_baseline_rmse_log10"].mean()),
                "train_term_prediction_rmse_log10": float(period_frame["train_term_prediction_rmse_log10"].mean()),
                "rmse_reduction_vs_zero_pct": float(period_frame["rmse_reduction_vs_zero_pct"].mean()),
            }
        )
        print(f"event holdout: period={spec.period_s:g}s", flush=True)
    return pd.DataFrame(rows)


def make_spatial_folds(frame: pd.DataFrame, n_splits: int, seed: int) -> np.ndarray:
    coordinates = frame[["lon", "lat"]].to_numpy(float)
    coordinates = np.nan_to_num(coordinates, nan=np.nanmedian(coordinates, axis=0))
    coordinates = StandardScaler().fit_transform(coordinates)
    return KMeans(n_clusters=n_splits, random_state=seed, n_init=30).fit_predict(coordinates)


def build_model(include_space: bool, seed: int) -> Pipeline:
    numeric_features = PHYSICAL_FEATURES + (SPATIAL_FEATURES if include_space else [])
    preprocess = ColumnTransformer(
        [
            (
                "numeric",
                SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
                numeric_features,
            ),
            (
                "categorical",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                    ]
                ),
                CATEGORICAL_FEATURES,
            ),
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


def cross_validate_station_models(
    station_terms: pd.DataFrame,
    min_records: int,
    n_splits: int,
    seed: int,
    target_model: str = "mf2013_site",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = station_terms[
        station_terms["model"].eq(target_model)
        & station_terms["n_records"].ge(min_records)
        & station_terms["lon"].notna()
        & station_terms["lat"].notna()
    ].copy()
    prediction_rows = []
    metric_rows = []
    for spec in PERIODS:
        frame = primary[primary["period_s"].eq(spec.period_s)].copy().reset_index(drop=True)
        y = frame["station_effect_log10"].to_numpy(float)
        weights = frame["n_records"].to_numpy(float)
        folds = make_spatial_folds(frame, n_splits=n_splits, seed=seed)
        for model_name, include_space in [("physical_hgb", False), ("physical_spatial_hgb", True)]:
            oof = np.full(len(frame), np.nan, dtype=float)
            for fold in sorted(np.unique(folds)):
                train = folds != fold
                test = folds == fold
                model = build_model(include_space=include_space, seed=seed + int(fold))
                model.fit(frame.loc[train], y[train], model__sample_weight=weights[train])
                oof[test] = model.predict(frame.loc[test])
            prediction_shift = float(np.average(oof, weights=weights))
            centered_oof = oof - prediction_shift
            for fold in sorted(np.unique(folds)):
                test = folds == fold
                fold_rmse = weighted_rmse(y[test], centered_oof[test], weights[test])
                fold_mae = weighted_mae(y[test], centered_oof[test], weights[test])
                fold_baseline_rmse = weighted_rmse(y[test], np.zeros(test.sum()), weights[test])
                fold_baseline_mae = weighted_mae(y[test], np.zeros(test.sum()), weights[test])
                metric_rows.append(
                    {
                        "period_s": spec.period_s,
                        "period_code": spec.code,
                        "residual_model": target_model,
                        "model": model_name,
                        "scope": "fold",
                        "fold": int(fold),
                        "n_stations": int(test.sum()),
                        "n_records_weight": int(weights[test].sum()),
                        "weighted_rmse_log10": fold_rmse,
                        "weighted_mae_log10": fold_mae,
                        "zero_baseline_rmse_log10": fold_baseline_rmse,
                        "zero_baseline_mae_log10": fold_baseline_mae,
                        "rmse_reduction_vs_zero_pct": 100.0 * (fold_baseline_rmse - fold_rmse) / fold_baseline_rmse,
                        "mae_reduction_vs_zero_pct": 100.0 * (fold_baseline_mae - fold_mae) / fold_baseline_mae,
                        "prediction_centering_shift_log10": prediction_shift,
                        "observed_predicted_correlation": float(np.corrcoef(y[test], centered_oof[test])[0, 1]),
                    }
                )
            zero = np.zeros_like(y)
            rmse = weighted_rmse(y, centered_oof, weights)
            mae = weighted_mae(y, centered_oof, weights)
            baseline_rmse = weighted_rmse(y, zero, weights)
            baseline_mae = weighted_mae(y, zero, weights)
            metric_rows.append(
                {
                    "period_s": spec.period_s,
                    "period_code": spec.code,
                    "residual_model": target_model,
                    "model": model_name,
                    "scope": "overall",
                    "fold": -1,
                    "n_stations": len(frame),
                    "n_records_weight": int(weights.sum()),
                    "weighted_rmse_log10": rmse,
                    "weighted_mae_log10": mae,
                    "zero_baseline_rmse_log10": baseline_rmse,
                    "zero_baseline_mae_log10": baseline_mae,
                    "rmse_reduction_vs_zero_pct": 100.0 * (baseline_rmse - rmse) / baseline_rmse,
                    "mae_reduction_vs_zero_pct": 100.0 * (baseline_mae - mae) / baseline_mae,
                    "prediction_centering_shift_log10": prediction_shift,
                    "observed_predicted_correlation": float(np.corrcoef(y, centered_oof)[0, 1]),
                }
            )
            block = frame[
                [
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
                ]
            ].copy()
            block["period_s"] = spec.period_s
            block["period_code"] = spec.code
            block["residual_model"] = target_model
            block["fold"] = folds
            block["model"] = model_name
            block["oof_prediction_log10"] = oof
            block["centered_oof_prediction_log10"] = centered_oof
            block["oof_multiplier"] = 10.0**centered_oof
            prediction_rows.append(block)
        print(f"models: period={spec.period_s:g}s stations={len(frame):,}", flush=True)
    return pd.concat(prediction_rows, ignore_index=True), pd.DataFrame(metric_rows)


def find_period_entry(archive: zipfile.ZipFile, period_code: str) -> str:
    candidates = []
    for name in archive.namelist():
        basename = Path(name).name.upper()
        if basename.endswith(".CSV") and period_code in basename and "PSV" not in basename:
            candidates.append(name)
    if not candidates:
        raise FileNotFoundError(f"No acceleration response map found for {period_code}")
    return sorted(candidates, key=lambda item: (len(item), item))[0]


def read_response_map_period(
    archive: zipfile.ZipFile,
    spec: PeriodSpec,
    meshcodes: set[str],
) -> pd.DataFrame:
    entry = find_period_entry(archive, spec.code)
    columns = ["CODE", *PROBABILITY_COLUMNS]
    selected = []
    with archive.open(entry) as handle:
        reader = pd.read_csv(handle, comment="#", header=None, names=columns, dtype=str, chunksize=500_000)
        for chunk in reader:
            chunk = chunk[chunk["CODE"].notna()].copy()
            chunk = chunk[~chunk["CODE"].str.upper().eq("CODE")]
            chunk["meshcode250_str"] = chunk["CODE"].map(normalize_mesh)
            chunk = chunk[chunk["meshcode250_str"].isin(meshcodes)]
            if chunk.empty:
                continue
            for col in PROBABILITY_COLUMNS:
                chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
            selected.append(chunk[["meshcode250_str", *PROBABILITY_COLUMNS]])
    out = pd.concat(selected, ignore_index=True)
    out = out.drop_duplicates("meshcode250_str", keep="first")
    out["period_s"] = spec.period_s
    out["period_code"] = spec.code
    return out


def build_hazard_values(
    predictions: pd.DataFrame,
    coefficients: pd.DataFrame,
    response_map: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    preferred = predictions[predictions["model"].eq("physical_spatial_hgb")].copy()
    preferred = preferred[preferred["meshcode250_str"].ne("")]
    meshcodes = set(preferred["meshcode250_str"])
    with zipfile.ZipFile(response_map) as archive:
        response = pd.concat(
            [read_response_map_period(archive, spec, meshcodes) for spec in PERIODS],
            ignore_index=True,
        )
    values = preferred.merge(
        response,
        on=["meshcode250_str", "period_s", "period_code"],
        how="inner",
        validate="many_to_one",
    )
    values = values.melt(
        id_vars=[column for column in values.columns if column not in PROBABILITY_COLUMNS],
        value_vars=list(PROBABILITY_COLUMNS),
        var_name="probability_column",
        value_name="official_vs400_sa_cm_s2",
    )
    values["probability_level"] = values["probability_column"].map(PROBABILITY_COLUMNS)
    values["official_vs400_sa_g"] = values["official_vs400_sa_cm_s2"] / G_IN_CM_S2
    ps_by_period = {spec.period_s: float(coefficients.loc[spec.coeff_key, "ps"]) for spec in PERIODS}
    vsmax_by_period = {spec.period_s: float(coefficients.loc[spec.coeff_key, "Vsmax"]) for spec in PERIODS}
    values["ps"] = values["period_s"].map(ps_by_period)
    values["vsmax"] = values["period_s"].map(vsmax_by_period)
    avs30 = pd.to_numeric(values["avs30"], errors="coerce")
    values["gs_station_log10"] = values["ps"] * np.log10(np.minimum(values["vsmax"], avs30) / 350.0)
    values["gs_vs400_log10"] = values["ps"] * np.log10(np.minimum(values["vsmax"], 400.0) / 350.0)
    values["surface_reference_factor"] = 10.0 ** (values["gs_station_log10"] - values["gs_vs400_log10"])
    values["ergodic_surface_sa_g"] = values["official_vs400_sa_g"] * values["surface_reference_factor"]
    values["adjusted_surface_sa_g"] = values["ergodic_surface_sa_g"] * values["oof_multiplier"]
    values["nonergodic_delta_pct"] = 100.0 * (values["oof_multiplier"] - 1.0)

    summary_rows = []
    for (period_s, probability), sub in values.groupby(["period_s", "probability_level"]):
        row = {
            "period_s": period_s,
            "probability_level": probability,
            "n_station_values": len(sub),
        }
        for col in [
            "official_vs400_sa_g",
            "surface_reference_factor",
            "ergodic_surface_sa_g",
            "oof_multiplier",
            "adjusted_surface_sa_g",
            "nonergodic_delta_pct",
        ]:
            for quantile in [0.05, 0.50, 0.95]:
                row[f"{col}_q{int(quantile * 100):02d}"] = float(sub[col].quantile(quantile))
        summary_rows.append(row)
    return values, pd.DataFrame(summary_rows)


def haversine_km(lon1: float, lat1: float, lon2: np.ndarray, lat2: np.ndarray) -> np.ndarray:
    radius_km = 6371.0088
    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    dphi = phi2 - phi1
    dlambda = np.deg2rad(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
    return 2.0 * radius_km * np.arcsin(np.sqrt(a))


def build_city_cases(hazard_values: pd.DataFrame) -> pd.DataFrame:
    stations = hazard_values[
        hazard_values["period_s"].eq(3.0) & hazard_values["probability_level"].eq("50y_10pct")
    ][["siteid2", "site_code", "network_label", "lon", "lat"]].drop_duplicates("siteid2")
    selected = []
    for city, (city_lon, city_lat) in CITIES.items():
        distances = haversine_km(
            city_lon,
            city_lat,
            stations["lon"].to_numpy(float),
            stations["lat"].to_numpy(float),
        )
        row = stations.iloc[int(np.argmin(distances))].to_dict()
        row.update({"city": city, "city_lon": city_lon, "city_lat": city_lat, "distance_km": float(np.min(distances))})
        selected.append(row)
    selected_frame = pd.DataFrame(selected)
    cases = hazard_values[
        hazard_values["probability_level"].eq("50y_10pct")
        & hazard_values["siteid2"].isin(selected_frame["siteid2"])
    ].merge(
        selected_frame[["siteid2", "city", "city_lon", "city_lat", "distance_km"]],
        on="siteid2",
        how="inner",
        validate="many_to_one",
    )
    return cases.sort_values(["city", "period_s"]).reset_index(drop=True)


def set_figure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def format_period_axis(ax: plt.Axes) -> None:
    ticks = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5]
    ax.set_xscale("log")
    ax.set_xticks(ticks, labels=["0.1", "0.2", "0.3", "0.5", "1", "2", "3", "5"])
    ax.set_xlabel("Period (s)")


def save_overview_figure(
    site: pd.DataFrame,
    source: pd.DataFrame,
    decomposition: pd.DataFrame,
    event_holdout: pd.DataFrame,
) -> None:
    set_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.15), constrained_layout=True)

    ax = axes[0, 0]
    stations = site[site["siteid2"].notna()].drop_duplicates("siteid2")
    ax.scatter(stations["lon"], stations["lat"], s=3, color="0.70", linewidths=0, label="K-NET/KiK-net")
    source_colors = {1: "#2A6F97", 2: "#D1495B", 3: "#E9C46A"}
    source_labels = {1: "Crustal", 2: "Interplate", 3: "Intraplate"}
    for source_type in [1, 2, 3]:
        sub = source[source["eq_location_type_id"].eq(source_type)]
        ax.scatter(
            sub["jem_lon"],
            sub["jem_lat"],
            s=8,
            color=source_colors[source_type],
            alpha=0.55,
            linewidths=0,
            label=source_labels[source_type],
        )
    ax.set_xlim(122, 151)
    ax.set_ylim(23, 47)
    ax.set_xlabel("Longitude (degrees E)")
    ax.set_ylabel("Latitude (degrees N)")
    ax.set_title("a  Strong-motion stations and earthquake sources")
    ax.legend(frameon=False, ncol=2, loc="lower left")

    metrics = decomposition.drop_duplicates(["period_s", "model"]).copy()
    basic = metrics[metrics["model"].eq("mf2013_basic")].set_index("period_s")
    site_model = metrics[metrics["model"].eq("mf2013_site")].set_index("period_s")
    common = sorted(set(basic.index) & set(site_model.index))
    mae_reduction = 100.0 * (basic.loc[common, "raw_mae_log10"] - site_model.loc[common, "raw_mae_log10"]) / basic.loc[common, "raw_mae_log10"]
    ax = axes[0, 1]
    ax.plot(common, mae_reduction, marker="o", color="#2A6F97", lw=1.4)
    ax.axhline(0.0, color="0.3", lw=0.8)
    format_period_axis(ax)
    ax.set_ylabel("MAE reduction from D1400/AVS30 (%)")
    ax.set_title("b  Direct-period MF2013 site-term gain")

    ax = axes[1, 0]
    site_metrics = site_model.loc[common]
    for column, label, color in [
        ("event_effect_weighted_std_log10", "Between-event", "#D1495B"),
        ("station_effect_weighted_std_log10", "Site-to-site", "#2A6F97"),
        ("remainder_rmse_log10", "Within-site", "0.35"),
    ]:
        ax.plot(common, site_metrics[column], marker="o", lw=1.3, color=color, label=label)
    format_period_axis(ax)
    ax.set_ylabel("Residual standard deviation (log10)")
    ax.set_title("c  Residual decomposition")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    repeat = event_holdout[event_holdout["scope"].eq("fold_mean")].sort_values("period_s")
    ax.plot(
        repeat["period_s"],
        repeat["train_test_station_correlation"],
        marker="o",
        color="#2A6F97",
        lw=1.4,
        label="Station-term correlation",
    )
    ax2 = ax.twinx()
    ax2.plot(
        repeat["period_s"],
        repeat["rmse_reduction_vs_zero_pct"],
        marker="s",
        color="#D1495B",
        lw=1.2,
        label="Fifth-group RMSE reduction",
    )
    format_period_axis(ax)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Four-to-fifth-group correlation")
    ax2.set_ylabel("Fifth-group RMSE reduction (%)")
    ax.set_title("d  Repeatability across earthquake groups")
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, frameon=False, loc="lower left")

    fig.savefig(OUT_OVERVIEW_PDF, bbox_inches="tight")
    fig.savefig(OUT_OVERVIEW_PNG, dpi=400, bbox_inches="tight")
    plt.close(fig)


def binned_station_effect(frame: pd.DataFrame, feature: str, bins: int = 12) -> pd.DataFrame:
    clean = frame[pd.to_numeric(frame[feature], errors="coerce").gt(0)].copy()
    clean["log_feature"] = np.log10(clean[feature].astype(float))
    clean["bin"] = pd.qcut(clean["log_feature"], q=bins, duplicates="drop")
    return clean.groupby("bin", observed=True).agg(
        x=("log_feature", "median"),
        y=("station_effect_log10", "median"),
        y25=("station_effect_log10", lambda values: values.quantile(0.25)),
        y75=("station_effect_log10", lambda values: values.quantile(0.75)),
    ).reset_index(drop=True)


def save_site_structure_figure(station_terms: pd.DataFrame, decomposition: pd.DataFrame) -> None:
    set_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.0), constrained_layout=True)
    model_style = {
        "mf2013_basic": ("MF2013 basic", "0.35"),
        "mf2013_site": ("MF2013 + D1400/AVS30", "#2A6F97"),
        "mf2013_site_ai_sensitivity": ("+ anomalous intensity", "#D1495B"),
    }
    for ax, feature, panel_title in [
        (axes[0, 0], "d1400", "a  Association with D1400"),
        (axes[0, 1], "dbase", "b  Association with basement depth"),
    ]:
        corr = decomposition[decomposition["association_feature"].eq(feature)]
        for model_name, (label, color) in model_style.items():
            sub = corr[corr["model"].eq(model_name)].sort_values("period_s")
            ax.plot(sub["period_s"], sub["station_effect_spearman"], marker="o", lw=1.3, color=color, label=label)
        ax.axhline(0.0, color="0.3", lw=0.8)
        format_period_axis(ax)
        ax.set_ylabel("Spearman correlation")
        ax.set_title(panel_title)
    axes[0, 1].legend(frameon=False, loc="best")

    ax = axes[1, 0]
    sa3 = station_terms[station_terms["period_s"].eq(3.0) & station_terms["n_records"].ge(20)]
    for model_name, label, color in [
        ("mf2013_basic", "Basic", "0.35"),
        ("mf2013_site", "D1400/AVS30", "#2A6F97"),
    ]:
        binned = binned_station_effect(sa3[sa3["model"].eq(model_name)], "d1400")
        ax.plot(binned["x"], binned["y"], marker="o", lw=1.4, color=color, label=label)
        ax.fill_between(binned["x"], binned["y25"], binned["y75"], color=color, alpha=0.14)
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_xlabel("log10 D1400 (m)")
    ax.set_ylabel("SA(3.0 s) station term (log10)")
    ax.set_title("c  Event-adjusted depth trend")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    compare = sa3[sa3["model"].isin(["mf2013_site", "mf2013_site_ai_sensitivity"])].pivot(
        index="siteid2", columns="model", values="station_effect_log10"
    ).dropna()
    x = compare["mf2013_site"].to_numpy(float)
    y = compare["mf2013_site_ai_sensitivity"].to_numpy(float)
    limit = float(np.quantile(np.abs(np.concatenate([x, y])), 0.995))
    ax.hexbin(x, y, gridsize=42, mincnt=1, cmap="Blues", linewidths=0)
    ax.plot([-limit, limit], [-limit, limit], color="0.25", lw=1.0, ls="--")
    ax.text(0.96, 0.04, "Dashed: equality", transform=ax.transAxes, ha="right", color="0.25")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("Primary station term (log10)")
    ax.set_ylabel("Anomalous-intensity station term (log10)")
    ax.set_title(f"d  Anomalous-intensity sensitivity (r = {np.corrcoef(x, y)[0, 1]:.3f})")

    fig.savefig(OUT_STRUCTURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_STRUCTURE_PNG, dpi=400, bbox_inches="tight")
    plt.close(fig)


def save_city_figure(city_cases: pd.DataFrame, predictions: pd.DataFrame) -> None:
    set_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.05), constrained_layout=True)
    colors = {
        "Tokyo": "#2A6F97",
        "Osaka": "#AA3377",
        "Sapporo": "#009988",
        "Sendai": "#EE7733",
        "Fukuoka": "#7A5195",
    }

    ax = axes[0, 0]
    sa3 = predictions[predictions["model"].eq("physical_spatial_hgb") & predictions["period_s"].eq(3.0)]
    ax.scatter(sa3["lon"], sa3["lat"], c=sa3["centered_oof_prediction_log10"], cmap="RdBu_r", s=5, alpha=0.45, linewidths=0)
    selected = city_cases.drop_duplicates("city")
    for _, row in selected.iterrows():
        ax.scatter(row["lon"], row["lat"], s=30, color=colors[row["city"]], edgecolor="black", linewidth=0.5, zorder=3)
        ax.text(row["lon"] + 0.25, row["lat"] + 0.12, row["city"], fontsize=7)
    ax.set_xlim(128, 146)
    ax.set_ylim(30, 46)
    ax.set_xlabel("Longitude (degrees E)")
    ax.set_ylabel("Latitude (degrees N)")
    ax.set_title("a  City-nearest matched stations")

    ax = axes[0, 1]
    for city in CITIES:
        sub = city_cases[city_cases["city"].eq(city)].sort_values("period_s")
        ax.plot(sub["period_s"], sub["oof_multiplier"], marker="o", lw=1.2, color=colors[city], label=city)
    ax.axhline(1.0, color="0.3", lw=0.8)
    format_period_axis(ax)
    ax.set_ylabel("Cross-validated station multiplier")
    ax.set_title("b  Direct-period station adjustments")
    ax.legend(frameon=False, ncol=2)

    ax = axes[1, 0]
    for city in CITIES:
        sub = city_cases[city_cases["city"].eq(city)].sort_values("period_s")
        ax.plot(sub["period_s"], sub["ergodic_surface_sa_g"], color=colors[city], lw=1.0, ls="--")
        ax.plot(sub["period_s"], sub["adjusted_surface_sa_g"], color=colors[city], lw=1.5, label=city)
    format_period_axis(ax)
    ax.set_yscale("log")
    ax.set_ylabel("50-year 10% spectral acceleration (g)")
    ax.set_title("c  Surface-reference response spectra")

    ax = axes[1, 1]
    long_period = city_cases[city_cases["period_s"].isin([1.0, 2.0, 3.0, 5.0])].copy()
    pivot = long_period.pivot(index="city", columns="period_s", values="nonergodic_delta_pct").reindex(CITIES)
    x = np.arange(len(pivot))
    width = 0.18
    period_colors = {1.0: "#A8DADC", 2.0: "#457B9D", 3.0: "#E76F51", 5.0: "#7A5195"}
    for index, period_s in enumerate([1.0, 2.0, 3.0, 5.0]):
        ax.bar(x + (index - 1.5) * width, pivot[period_s], width=width, color=period_colors[period_s], label=f"{period_s:g} s")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_xticks(x, pivot.index, rotation=20, ha="right")
    ax.set_ylabel("Spectral ordinate change (%)")
    ax.set_title("d  Long-period changes")
    ax.legend(frameon=False, ncol=2)

    fig.savefig(OUT_CITY_PDF, bbox_inches="tight")
    fig.savefig(OUT_CITY_PNG, dpi=400, bbox_inches="tight")
    plt.close(fig)


def save_figure(
    station_terms: pd.DataFrame,
    model_metrics: pd.DataFrame,
    predictions: pd.DataFrame,
    hazard_summary: pd.DataFrame,
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    set_figure_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.2), constrained_layout=True)

    ax = axes[0, 0]
    map_data = station_terms[
        station_terms["model"].eq("mf2013_site")
        & station_terms["period_s"].eq(3.0)
        & station_terms["n_records"].ge(20)
    ].copy()
    vmax = float(np.nanquantile(np.abs(map_data["station_effect_log10"]), 0.98))
    scatter = ax.scatter(
        map_data["lon"],
        map_data["lat"],
        c=map_data["station_effect_log10"],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        s=7,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xlim(128, 146)
    ax.set_ylim(30, 46)
    ax.set_xlabel("Longitude (degrees E)")
    ax.set_ylabel("Latitude (degrees N)")
    ax.set_title("a  Event-adjusted SA(3.0 s) station terms")
    fig.colorbar(scatter, ax=ax, label="Station term (log10)", fraction=0.046, pad=0.03)

    ax = axes[0, 1]
    colors = {"physical_hgb": "#2A6F97", "physical_spatial_hgb": "#D1495B"}
    labels = {"physical_hgb": "Site parameters", "physical_spatial_hgb": "Site + regional"}
    for model_name in ["physical_hgb", "physical_spatial_hgb"]:
        sub = model_metrics[
            model_metrics["model"].eq(model_name) & model_metrics["scope"].eq("overall")
        ].sort_values("period_s")
        ax.plot(
            sub["period_s"],
            sub["rmse_reduction_vs_zero_pct"],
            marker="o",
            lw=1.4,
            color=colors[model_name],
            label=labels[model_name],
        )
    ax.axhline(0.0, color="0.3", lw=0.8)
    format_period_axis(ax)
    ax.set_ylabel("Spatial-block RMSE reduction (%)")
    ax.set_title("b  Spatial estimation of station terms")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    pred = predictions[predictions["model"].eq("physical_spatial_hgb")].copy()
    rows = []
    for period_s, sub in pred.groupby("period_s"):
        rows.append(
            {
                "period_s": period_s,
                "q05": sub["oof_multiplier"].quantile(0.05),
                "q50": sub["oof_multiplier"].quantile(0.50),
                "q95": sub["oof_multiplier"].quantile(0.95),
            }
        )
    quantiles = pd.DataFrame(rows).sort_values("period_s")
    ax.plot(quantiles["period_s"], quantiles["q50"], color="#2A6F97", marker="o", lw=1.4)
    ax.fill_between(quantiles["period_s"], quantiles["q05"], quantiles["q95"], color="#2A6F97", alpha=0.2)
    ax.axhline(1.0, color="0.3", lw=0.8)
    format_period_axis(ax)
    ax.set_ylabel("Cross-validated multiplier")
    ax.set_title("c  Zero-centred station adjustment")

    ax = axes[1, 1]
    spectrum = hazard_summary[hazard_summary["probability_level"].eq("50y_10pct")].sort_values("period_s")
    ax.plot(spectrum["period_s"], spectrum["official_vs400_sa_g_q50"], color="0.25", marker="o", lw=1.2, label="Official Vs400")
    ax.plot(spectrum["period_s"], spectrum["ergodic_surface_sa_g_q50"], color="#2A6F97", marker="o", lw=1.2, label="MF2013 surface")
    ax.plot(spectrum["period_s"], spectrum["adjusted_surface_sa_g_q50"], color="#D1495B", marker="o", lw=1.2, label="Station adjusted")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5], labels=["0.1", "0.2", "0.3", "0.5", "1", "2", "3", "5"])
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Median spectral acceleration (g)")
    ax.set_title("d  Matched-station 50-year 10% spectra")
    ax.legend(frameon=False)

    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=400, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    station_terms: pd.DataFrame,
    decomposition: pd.DataFrame,
    model_metrics: pd.DataFrame,
    predictions: pd.DataFrame,
    event_holdout: pd.DataFrame,
    hazard_summary: pd.DataFrame,
) -> None:
    primary_metrics = model_metrics[
        model_metrics["model"].eq("physical_spatial_hgb") & model_metrics["scope"].eq("overall")
    ].sort_values("period_s")
    sa3 = primary_metrics[primary_metrics["period_s"].eq(3.0)].iloc[0]
    pred3 = predictions[
        predictions["model"].eq("physical_spatial_hgb") & predictions["period_s"].eq(3.0)
    ]
    hz3 = hazard_summary[
        hazard_summary["period_s"].eq(3.0) & hazard_summary["probability_level"].eq("50y_10pct")
    ].iloc[0]
    ai_compare = station_terms[
        station_terms["period_s"].eq(3.0)
        & station_terms["model"].isin(["mf2013_site", "mf2013_site_ai_sensitivity"])
        & station_terms["n_records"].ge(20)
    ].pivot(index="siteid2", columns="model", values="station_effect_log10").dropna()
    ai_corr = float(ai_compare.corr().iloc[0, 1])
    repeat3 = event_holdout[event_holdout["period_s"].eq(3.0) & event_holdout["scope"].eq("fold_mean")].iloc[0]

    lines = [
        "# Event-adjusted direct-period station analysis",
        "",
        "## Primary design",
        "",
        "- RotD100 response spectra at eight periods (0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, and 5.0 s) are read directly from the public flatfile to match the MF2013 horizontal-vector response definition.",
        "- Each residual is decomposed into a global intercept, a zero-centred event term, a zero-centred station term, and a record remainder.",
        "- Station prediction uses five spatial blocks and only out-of-fold predictions in the hazard sensitivity calculation.",
        "- Official response-map ordinates are retained at Vs=400 m/s and are separately projected to each station AVS30 before the station adjustment is applied.",
        "",
        "## SA(3.0 s) result",
        "",
        f"- Supported stations: {int(sa3['n_stations']):,}.",
        f"- Spatial-block RMSE reduction versus a zero station term: {sa3['rmse_reduction_vs_zero_pct']:.1f}%.",
        f"- Observed-predicted station-term correlation: {sa3['observed_predicted_correlation']:.3f}.",
        f"- Across held-out event groups, the mean train-test station-term correlation is {repeat3['train_test_station_correlation']:.3f} and the mean test RMSE reduction is {repeat3['rmse_reduction_vs_zero_pct']:.1f}%.",
        f"- Cross-validated station multiplier q05/q50/q95: {pred3['oof_multiplier'].quantile(0.05):.3f} / {pred3['oof_multiplier'].quantile(0.50):.3f} / {pred3['oof_multiplier'].quantile(0.95):.3f}.",
        f"- Matched-station 50-year 10% official Vs400 median: {hz3['official_vs400_sa_g_q50']:.3f} g.",
        f"- MF2013 surface-reference median: {hz3['ergodic_surface_sa_g_q50']:.3f} g.",
        f"- Station-adjusted surface-reference median: {hz3['adjusted_surface_sa_g_q50']:.3f} g.",
        "",
        "## Model boundary",
        "",
        "- The primary attenuation backbone contains the MF2013 basic, D1400, and AVS30 terms. It is labelled as such and is not described as the complete official implementation.",
        f"- The rule-based anomalous-intensity sensitivity gives an SA(3.0 s) station-term correlation of {ai_corr:.3f} with the primary decomposition.",
        "- PH is an event-constant period term for qualifying Philippine Sea Plate intraplate earthquakes. The event fixed effect absorbs it for station-term estimation; the global intercept is not interpreted as an official MF2013 bias.",
        "- The response-map calculation is a matched-station surface-spectrum sensitivity analysis, not an official source-level J-SHIS PSHA rerun.",
        "",
        "## Output files",
        "",
        f"- `{OUT_STATION_TERMS.relative_to(ROOT)}`",
        f"- `{OUT_DECOMPOSITION.relative_to(ROOT)}`",
        f"- `{OUT_MODEL_PREDICTIONS.relative_to(ROOT)}`",
        f"- `{OUT_MODEL_METRICS.relative_to(ROOT)}`",
        f"- `{OUT_EVENT_HOLDOUT.relative_to(ROOT)}`",
        f"- `{OUT_HAZARD_VALUES.relative_to(ROOT)}`",
        f"- `{OUT_HAZARD_SUMMARY.relative_to(ROOT)}`",
        f"- `{OUT_CITY_CASES.relative_to(ROOT)}`",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    coefficients = load_coefficients(args.coefficients)
    site, source = load_metadata(args.flatfile)
    records = load_record_residuals(args.flatfile, coefficients, site, source, args.chunksize)
    station_terms, decomposition = decompose_all_periods(records, site)
    audit_response_components(records, station_terms, decomposition)
    event_holdout = event_holdout_repeatability(
        records,
        n_splits=args.spatial_folds,
        seed=args.seed,
    )
    predictions, model_metrics = cross_validate_station_models(
        station_terms,
        min_records=args.min_station_records,
        n_splits=args.spatial_folds,
        seed=args.seed,
    )
    station_terms.to_csv(OUT_STATION_TERMS, index=False)
    decomposition.to_csv(OUT_DECOMPOSITION, index=False)
    predictions.to_csv(OUT_MODEL_PREDICTIONS, index=False)
    model_metrics.to_csv(OUT_MODEL_METRICS, index=False)
    event_holdout.to_csv(OUT_EVENT_HOLDOUT, index=False)
    hazard_values, hazard_summary = build_hazard_values(predictions, coefficients, args.response_map)
    city_cases = build_city_cases(hazard_values)
    hazard_values.to_csv(OUT_HAZARD_VALUES, index=False)
    hazard_summary.to_csv(OUT_HAZARD_SUMMARY, index=False)
    city_cases.to_csv(OUT_CITY_CASES, index=False)
    save_overview_figure(site, source, decomposition, event_holdout)
    save_site_structure_figure(station_terms, decomposition)
    save_figure(station_terms, model_metrics, predictions, hazard_summary)
    save_city_figure(city_cases, predictions)
    write_audit(station_terms, decomposition, model_metrics, predictions, event_holdout, hazard_summary)
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=DEFAULT_FLATFILE)
    parser.add_argument("--coefficients", type=Path, default=DEFAULT_COEFFICIENTS)
    parser.add_argument("--response-map", type=Path, default=DEFAULT_RESPONSE_MAP)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--min-station-records", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
