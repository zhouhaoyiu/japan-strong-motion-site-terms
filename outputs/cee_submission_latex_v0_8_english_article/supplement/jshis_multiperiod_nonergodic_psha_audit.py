#!/usr/bin/env python3
"""All-period official-response-spectrum non-ergodic station-term audit.

This script extends the previous SA(3.0 s)-only hazard-facing check to the
available J-SHIS response-spectrum periods. It builds station and sampled-grid
correction terms at the periods supported by the residual model, interpolates
between those anchors, and applies the resulting multipliers to official
J-SHIS 250 m response-spectrum map ordinates.

The calculation is a production-style official-map propagation audit. It is
not an official source-level J-SHIS PSHA rerun because the local package does
not contain the complete official source recurrence, period-dependent sigma,
spatial correlation, and logic-tree implementation.
"""

from __future__ import annotations

import argparse
import math
import shutil
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
ARTICLE_DIR = OUTPUT_DIR / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT_DIR = ARTICLE_DIR / "supplement"
FIGURE_DIR = ARTICLE_DIR / "figures"

PREDICTIONS = OUTPUT_DIR / "jshis_station_model_ablation_predictions.csv"
SITE_SCHEMA = OUTPUT_DIR / "jshis_site_schema_v2024_sub1.csv"
MESH_LINK = SUPPLEMENT_DIR / "jshis_public_psha_10467_nearest_station_borehole_mesh_link.csv"
RESP_ZIP = ROOT / "work" / "external_data" / "jshis_respmap" / "P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip"

OUT_STATION_CORR = OUTPUT_DIR / "jshis_multiperiod_station_corrections.csv"
OUT_VALIDATION = OUTPUT_DIR / "jshis_multiperiod_station_surface_validation.csv"
OUT_VALIDATION_SUMMARY = OUTPUT_DIR / "jshis_multiperiod_station_surface_validation_summary.csv"
OUT_MESH = OUTPUT_DIR / "jshis_multiperiod_continuous_correction_surface_mesh.csv"
OUT_STATION_VALUES = OUTPUT_DIR / "jshis_multiperiod_official_response_station_values.csv"
OUT_MESH_VALUES = OUTPUT_DIR / "jshis_multiperiod_official_response_mesh_values.csv"
OUT_SUMMARY = OUTPUT_DIR / "jshis_multiperiod_nonergodic_psha_summary.csv"
OUT_AUDIT = OUTPUT_DIR / "jshis_multiperiod_nonergodic_psha_audit.md"
OUT_FIGURE = FIGURE_DIR / "figure17_multiperiod_nonergodic_psha.png"

G_IN_CM_S2 = 980.665
K_NEIGHBORS = 24
IDW_POWER = 2.0
EPS_KM = 1.0e-3

PERIOD_CODES = {
    "P0010": 0.1,
    "P0020": 0.2,
    "P0030": 0.3,
    "P0050": 0.5,
    "P0100": 1.0,
    "P0200": 2.0,
    "P0300": 3.0,
    "P0500": 5.0,
}
PROBABILITY_COLUMNS = {
    "T50_P02_SA": "50y_2pct",
    "T50_P05_SA": "50y_5pct",
    "T50_P10_SA": "50y_10pct",
    "T50_P39_SA": "50y_39pct",
}
ANCHOR_TARGETS = {
    "SA(0.3s) RotD50": {"key": "sa03", "period_s": 0.3},
    "SA(1.0s) RotD50": {"key": "sa10", "period_s": 1.0},
    "SA(3.0s) RotD50": {"key": "sa30", "period_s": 3.0},
}
OPTIONAL_TARGETS = {
    "PGA RotD50": {"key": "pga", "period_s": np.nan},
    **ANCHOR_TARGETS,
}

MODEL_FEATURES = [
    "lon",
    "lat",
    "elevation",
    "sensor_depth_glminus",
    "dist_vf_mf13_nejapan",
    "dist_vf_mf13_swjapan",
    "log_vs10",
    "log_vs20",
    "log_vs30",
    "log_avs30",
    "log_d1100",
    "log_d1400",
    "log_d1700",
    "log_d2100",
    "log_dbase",
]
SITE_SCHEMA_COLUMNS = [
    "siteid2",
    "site_code",
    "lon",
    "lat",
    "obs_network_id",
    "installation_situation_id",
    "elevation",
    "sensor_depth_glminus",
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


def normalize_mesh(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def lonlat_to_xy_km(lon: pd.Series | np.ndarray, lat: pd.Series | np.ndarray, ref_lat: float | None = None) -> np.ndarray:
    lon_arr = np.asarray(lon, dtype=float)
    lat_arr = np.asarray(lat, dtype=float)
    lat0 = float(np.nanmedian(lat_arr) if ref_lat is None else ref_lat)
    x = lon_arr * 111.32 * np.cos(np.deg2rad(lat0))
    y = lat_arr * 110.57
    return np.column_stack([x, y])


def weighted_mae(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(np.abs(y_true - y_pred), weights=weights))


def weighted_rmse(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.average((y_true - y_pred) ** 2, weights=weights)))


def idw_predict(
    train_xy: np.ndarray,
    train_y: np.ndarray,
    query_xy: np.ndarray,
    k_neighbors: int = K_NEIGHBORS,
    power: float = IDW_POWER,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    k = min(k_neighbors, len(train_y))
    nn = NearestNeighbors(n_neighbors=k, algorithm="ball_tree")
    nn.fit(train_xy)
    dist, idx = nn.kneighbors(query_xy, return_distance=True)
    close = dist <= EPS_KM
    weights = 1.0 / np.maximum(dist, EPS_KM) ** power
    pred = np.empty(len(query_xy), dtype=float)
    local_std = np.empty(len(query_xy), dtype=float)
    for i in range(len(query_xy)):
        if close[i].any():
            vals = train_y[idx[i, close[i]]]
            pred[i] = float(np.mean(vals))
            local_std[i] = 0.0 if len(vals) == 1 else float(np.std(vals))
            continue
        vals = train_y[idx[i]]
        w = weights[i]
        pred[i] = float(np.average(vals, weights=w))
        local_std[i] = float(np.sqrt(np.average((vals - pred[i]) ** 2, weights=w)))
    return pred, dist[:, 0], local_std


def add_log_features(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for col in ["vs10", "vs20", "vs30", "avs30", "d1100", "d1400", "d1700", "d2100", "dbase"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
        out[f"log_{col}"] = np.log10(out[col].clip(lower=1.0))
    return out


def make_site_space_model(seed: int) -> Pipeline:
    numeric = Pipeline([("imputer", SimpleImputer(strategy="median"))])
    return Pipeline(
        [
            ("preprocess", ColumnTransformer([("numeric", numeric, MODEL_FEATURES)])),
            (
                "model",
                HistGradientBoostingRegressor(
                    loss="squared_error",
                    learning_rate=0.04,
                    max_iter=400,
                    max_leaf_nodes=15,
                    min_samples_leaf=35,
                    l2_regularization=0.1,
                    random_state=seed,
                ),
            ),
        ]
    )


def load_station_corrections() -> pd.DataFrame:
    pred = pd.read_csv(PREDICTIONS)
    pred = pred[pred["split"].eq("spatial_block") & pred["target_label"].isin(OPTIONAL_TARGETS)].copy()
    pred["target_key"] = pred["target_label"].map(lambda x: OPTIONAL_TARGETS[x]["key"])
    pred["anchor_period_s"] = pred["target_label"].map(lambda x: OPTIONAL_TARGETS[x]["period_s"])
    pred = pred.rename(columns={"full_site_location_prediction": "correction_log10"})

    schema = pd.read_csv(SITE_SCHEMA, usecols=SITE_SCHEMA_COLUMNS)
    schema = schema.sort_values(["siteid2", "site_code"]).drop_duplicates("siteid2", keep="last")
    stations = pred.merge(schema, on=["siteid2", "site_code"], how="left", validate="many_to_one")
    numeric = [
        "lon",
        "lat",
        "n_records",
        "correction_log10",
        "fold",
        "elevation",
        "sensor_depth_glminus",
        "dist_vf_mf13_nejapan",
        "dist_vf_mf13_swjapan",
        "vs10",
        "vs20",
        "vs30",
        "avs30",
        "d1100",
        "d1400",
        "d1700",
        "d2100",
        "dbase",
    ]
    for col in numeric:
        stations[col] = pd.to_numeric(stations[col], errors="coerce")
    stations = stations.dropna(subset=["lon", "lat", "n_records", "correction_log10", "fold"]).copy()
    stations["fold"] = stations["fold"].astype(int)
    stations["n_records"] = stations["n_records"].clip(lower=1)
    stations["meshcode250_str"] = stations["meshcode250"].map(normalize_mesh)
    stations["multiplier"] = 10.0 ** stations["correction_log10"].astype(float)
    return add_log_features(stations)


def validate_surfaces(stations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref_lat = float(stations["lat"].median())
    rows = []
    summary_rows = []
    for target_label, sub in stations[stations["target_label"].isin(ANCHOR_TARGETS)].groupby("target_label", sort=False):
        sub = sub.reset_index(drop=True)
        xy = lonlat_to_xy_km(sub["lon"], sub["lat"], ref_lat=ref_lat)
        y = sub["correction_log10"].to_numpy(float)
        weights = sub["n_records"].to_numpy(float)
        for fold in sorted(sub["fold"].unique()):
            test = sub["fold"].eq(fold).to_numpy()
            train = ~test
            idw_pred, d1, local_std = idw_predict(xy[train], y[train], xy[test])
            model = make_site_space_model(seed=20260617 + int(fold))
            model.fit(sub.loc[train, MODEL_FEATURES], y[train], model__sample_weight=weights[train])
            model_pred = model.predict(sub.loc[test, MODEL_FEATURES])
            train_mean = float(np.average(y[train], weights=weights[train]))
            baseline_zero = np.zeros_like(y[test])
            baseline_mean = np.full_like(y[test], train_mean)
            for method, pred_values in [("idw", idw_pred), ("site_space_hgb", model_pred)]:
                mae = weighted_mae(y[test], pred_values, weights[test])
                rmse = weighted_rmse(y[test], pred_values, weights[test])
                zero_mae = weighted_mae(y[test], baseline_zero, weights[test])
                zero_rmse = weighted_rmse(y[test], baseline_zero, weights[test])
                mean_mae = weighted_mae(y[test], baseline_mean, weights[test])
                mean_rmse = weighted_rmse(y[test], baseline_mean, weights[test])
                summary_rows.append(
                    {
                        "target_label": target_label,
                        "target_key": ANCHOR_TARGETS[target_label]["key"],
                        "method": method,
                        "fold": int(fold),
                        "n_stations": int(test.sum()),
                        "weighted_mae": mae,
                        "weighted_rmse": rmse,
                        "zero_baseline_mae": zero_mae,
                        "zero_baseline_rmse": zero_rmse,
                        "mean_baseline_mae": mean_mae,
                        "mean_baseline_rmse": mean_rmse,
                        "rmse_reduction_vs_zero_pct": (zero_rmse - rmse) / zero_rmse * 100.0 if zero_rmse else np.nan,
                        "rmse_reduction_vs_mean_pct": (mean_rmse - rmse) / mean_rmse * 100.0 if mean_rmse else np.nan,
                    }
                )
            block = sub.loc[test, ["siteid2", "site_code", "network_label", "target_label", "target_key", "fold", "n_records", "lon", "lat"]].copy()
            block["observed_correction_log10"] = y[test]
            block["idw_predicted_correction_log10"] = idw_pred
            block["site_space_predicted_correction_log10"] = model_pred
            block["nearest_train_distance_km"] = d1
            block["local_idw_std_log10"] = local_std
            rows.append(block)

        pred_frame = pd.concat([r for r in rows if r["target_label"].iloc[0] == target_label], ignore_index=True)
        for method, col in [("idw", "idw_predicted_correction_log10"), ("site_space_hgb", "site_space_predicted_correction_log10")]:
            all_y = pred_frame["observed_correction_log10"].to_numpy(float)
            all_pred = pred_frame[col].to_numpy(float)
            all_weights = pred_frame["n_records"].to_numpy(float)
            zero_rmse = weighted_rmse(all_y, np.zeros_like(all_y), all_weights)
            mean_rmse = weighted_rmse(all_y, np.full_like(all_y, np.average(all_y, weights=all_weights)), all_weights)
            rmse = weighted_rmse(all_y, all_pred, all_weights)
            summary_rows.append(
                {
                    "target_label": target_label,
                    "target_key": ANCHOR_TARGETS[target_label]["key"],
                    "method": f"overall_{method}",
                    "fold": -1,
                    "n_stations": int(len(pred_frame)),
                    "weighted_mae": weighted_mae(all_y, all_pred, all_weights),
                    "weighted_rmse": rmse,
                    "zero_baseline_mae": weighted_mae(all_y, np.zeros_like(all_y), all_weights),
                    "zero_baseline_rmse": zero_rmse,
                    "mean_baseline_mae": weighted_mae(all_y, np.full_like(all_y, np.average(all_y, weights=all_weights)), all_weights),
                    "mean_baseline_rmse": mean_rmse,
                    "rmse_reduction_vs_zero_pct": (zero_rmse - rmse) / zero_rmse * 100.0 if zero_rmse else np.nan,
                    "rmse_reduction_vs_mean_pct": (mean_rmse - rmse) / mean_rmse * 100.0 if mean_rmse else np.nan,
                }
            )
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summary_rows)


def build_mesh_features() -> pd.DataFrame:
    mesh = pd.read_csv(MESH_LINK)
    for col in ["lon", "lat"]:
        mesh[col] = pd.to_numeric(mesh[col], errors="coerce")
    mesh = mesh.dropna(subset=["lon", "lat"]).copy()
    mesh["meshcode250_str"] = mesh["meshcode250"].map(normalize_mesh)

    schema = pd.read_csv(SITE_SCHEMA, usecols=SITE_SCHEMA_COLUMNS)
    schema = schema.sort_values(["siteid2", "site_code"]).drop_duplicates("siteid2", keep="last")
    nearest = schema.drop(columns=["site_code", "lon", "lat", "meshcode250", "meshcode3"]).add_prefix("nearest_schema_")
    mesh = mesh.merge(nearest, left_on="nearest_siteid2", right_on="nearest_schema_siteid2", how="left", validate="many_to_one")
    for col in [
        "elevation",
        "sensor_depth_glminus",
        "dist_vf_mf13_nejapan",
        "dist_vf_mf13_swjapan",
        "vs10",
        "vs20",
        "vs30",
        "avs30",
        "d1100",
        "d1400",
        "d1700",
        "d2100",
        "dbase",
    ]:
        mesh[col] = pd.to_numeric(mesh.get(f"nearest_schema_{col}"), errors="coerce")
    return add_log_features(mesh)


def build_mesh_surface(stations: pd.DataFrame) -> pd.DataFrame:
    mesh = build_mesh_features()
    ref_lat = float(pd.concat([stations["lat"], mesh["lat"]]).median())
    mesh_xy = lonlat_to_xy_km(mesh["lon"], mesh["lat"], ref_lat=ref_lat)
    out = mesh.copy()
    for i, (target_label, info) in enumerate(ANCHOR_TARGETS.items()):
        key = info["key"]
        sub = stations[stations["target_label"].eq(target_label)].copy()
        sub_xy = lonlat_to_xy_km(sub["lon"], sub["lat"], ref_lat=ref_lat)
        y = sub["correction_log10"].to_numpy(float)
        idw_pred, d1, local_std = idw_predict(sub_xy, y, mesh_xy)
        model = make_site_space_model(seed=20260617 + i)
        model.fit(sub[MODEL_FEATURES], y, model__sample_weight=sub["n_records"].to_numpy(float))
        model_pred = model.predict(out[MODEL_FEATURES])
        out[f"idw_{key}_correction_log10"] = idw_pred
        out[f"idw_{key}_multiplier"] = 10.0**idw_pred
        out[f"surface_{key}_correction_log10"] = model_pred
        out[f"surface_{key}_multiplier"] = 10.0**model_pred
        out[f"surface_{key}_nearest_station_distance_km"] = d1
        out[f"surface_{key}_local_idw_std_log10"] = local_std
    return out


def correction_source(period_s: float) -> str:
    if np.isclose(period_s, 0.3) or np.isclose(period_s, 1.0) or np.isclose(period_s, 3.0):
        return "direct_anchor"
    if period_s < 0.3:
        return "clamped_to_0.3s"
    if period_s > 3.0:
        return "clamped_to_3.0s"
    return "log_period_interpolated"


def interpolate_corrections(anchor_frame: pd.DataFrame, prefix: str = "") -> pd.DataFrame:
    rows = []
    anchors = np.array([0.3, 1.0, 3.0], dtype=float)
    log_anchors = np.log10(anchors)
    values = anchor_frame[[f"{prefix}sa03_correction_log10", f"{prefix}sa10_correction_log10", f"{prefix}sa30_correction_log10"]].to_numpy(float)
    for period_code, period_s in PERIOD_CODES.items():
        logp = math.log10(period_s)
        interp = np.array([np.interp(logp, log_anchors, row, left=row[0], right=row[-1]) for row in values], dtype=float)
        block = anchor_frame[["meshcode250_str"]].copy()
        for col in ["siteid2", "site_code", "network_label", "lon", "lat"]:
            if col in anchor_frame:
                block[col] = anchor_frame[col].values
        block["period_code"] = period_code
        block["period_s"] = period_s
        block["correction_log10"] = interp
        block["multiplier"] = 10.0**interp
        block["correction_source"] = correction_source(period_s)
        rows.append(block)
    return pd.concat(rows, ignore_index=True)


def build_station_anchor_wide(stations: pd.DataFrame) -> pd.DataFrame:
    anchor = stations[stations["target_label"].isin(ANCHOR_TARGETS)].copy()
    pieces = []
    base_cols = ["siteid2", "site_code", "network_label", "lon", "lat", "meshcode250_str", "n_records", "fold"]
    base = anchor.sort_values("target_key").drop_duplicates(["siteid2", "site_code"])[base_cols].copy()
    for target_label, info in ANCHOR_TARGETS.items():
        sub = anchor[anchor["target_label"].eq(target_label)][["siteid2", "site_code", "correction_log10", "multiplier"]].copy()
        sub = sub.rename(
            columns={
                "correction_log10": f"{info['key']}_correction_log10",
                "multiplier": f"{info['key']}_multiplier",
            }
        )
        pieces.append(sub)
    out = base.copy()
    for piece in pieces:
        out = out.merge(piece, on=["siteid2", "site_code"], how="left", validate="one_to_one")
    return out.dropna(subset=["sa03_correction_log10", "sa10_correction_log10", "sa30_correction_log10"])


def build_mesh_anchor_wide(mesh: pd.DataFrame) -> pd.DataFrame:
    out = mesh[
        [
            "meshcode250_str",
            "meshcode250",
            "meshcode3",
            "lon",
            "lat",
            "is_tail",
            "is_ne_hokkaido",
            "basin_period_proxy_1_3s",
            "nearest_station_distance_km",
            "nearest_siteid2",
            "nearest_site_code",
        ]
    ].copy()
    for key in ["sa03", "sa10", "sa30"]:
        out[f"{key}_correction_log10"] = mesh[f"surface_{key}_correction_log10"]
        out[f"{key}_multiplier"] = mesh[f"surface_{key}_multiplier"]
    return out


def find_period_entry(zf: zipfile.ZipFile, period_code: str) -> str:
    candidates = []
    for name in zf.namelist():
        upper = Path(name).name.upper()
        if upper.endswith(".CSV") and period_code.upper() in upper:
            candidates.append(name)
    if not candidates:
        raise FileNotFoundError(f"No response-spectrum CSV found for {period_code}")
    return sorted(candidates, key=lambda item: (("PSV" in item.upper()), len(item), item))[0]


def read_response_period(zf: zipfile.ZipFile, period_code: str, meshcodes: set[str]) -> pd.DataFrame:
    entry = find_period_entry(zf, period_code)
    names = ["CODE", "T50_P02_SA", "T50_P05_SA", "T50_P10_SA", "T50_P39_SA"]
    selected = []
    with zf.open(entry) as handle:
        chunks = pd.read_csv(handle, comment="#", header=None, names=names, dtype=str, chunksize=500_000)
        for chunk in chunks:
            chunk = chunk[chunk["CODE"].notna()].copy()
            chunk = chunk[~chunk["CODE"].astype(str).str.upper().eq("CODE")]
            chunk["meshcode250_str"] = chunk["CODE"].map(normalize_mesh)
            chunk = chunk[chunk["meshcode250_str"].isin(meshcodes)]
            if chunk.empty:
                continue
            for col in names[1:]:
                chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
            chunk["period_code"] = period_code
            chunk["period_s"] = PERIOD_CODES[period_code]
            selected.append(chunk[["meshcode250_str", "period_code", "period_s", *names[1:]]])
    if not selected:
        return pd.DataFrame(columns=["meshcode250_str", "period_code", "period_s", *names[1:]])
    return pd.concat(selected, ignore_index=True)


def load_response_maps(meshcodes: set[str]) -> pd.DataFrame:
    with zipfile.ZipFile(RESP_ZIP) as zf:
        frames = [read_response_period(zf, code, meshcodes) for code in PERIOD_CODES]
    response = pd.concat(frames, ignore_index=True)
    response = response.melt(
        id_vars=["meshcode250_str", "period_code", "period_s"],
        value_vars=list(PROBABILITY_COLUMNS),
        var_name="probability_column",
        value_name="official_sa_cm_s2",
    )
    response["probability_level"] = response["probability_column"].map(PROBABILITY_COLUMNS)
    response["official_sa_g"] = pd.to_numeric(response["official_sa_cm_s2"], errors="coerce") / G_IN_CM_S2
    return response.dropna(subset=["official_sa_g"]).copy()


def apply_response_maps(entity: pd.DataFrame, corrections: pd.DataFrame, scope: str) -> pd.DataFrame:
    meshcodes = set(entity["meshcode250_str"].dropna().astype(str))
    response = load_response_maps(meshcodes)
    merged = response.merge(corrections, on=["meshcode250_str", "period_code", "period_s"], how="inner", validate="many_to_many")
    merged["scope"] = scope
    merged["corrected_sa_g"] = merged["official_sa_g"] * merged["multiplier"]
    merged["delta_sa_g"] = merged["corrected_sa_g"] - merged["official_sa_g"]
    merged["delta_pct"] = 100.0 * (merged["multiplier"] - 1.0)
    return merged


def build_summary(station_values: pd.DataFrame, mesh_values: pd.DataFrame, validation_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scope, frame in [("station", station_values), ("sampled_grid", mesh_values)]:
        for (prob, period_s, source), sub in frame.groupby(["probability_level", "period_s", "correction_source"], sort=True):
            rows.append(
                {
                    "scope": scope,
                    "probability_level": prob,
                    "period_s": float(period_s),
                    "correction_source": source,
                    "n_values": int(len(sub)),
                    "official_sa_g_q05": float(sub["official_sa_g"].quantile(0.05)),
                    "official_sa_g_q50": float(sub["official_sa_g"].quantile(0.50)),
                    "official_sa_g_q95": float(sub["official_sa_g"].quantile(0.95)),
                    "corrected_sa_g_q05": float(sub["corrected_sa_g"].quantile(0.05)),
                    "corrected_sa_g_q50": float(sub["corrected_sa_g"].quantile(0.50)),
                    "corrected_sa_g_q95": float(sub["corrected_sa_g"].quantile(0.95)),
                    "multiplier_q05": float(sub["multiplier"].quantile(0.05)),
                    "multiplier_q50": float(sub["multiplier"].quantile(0.50)),
                    "multiplier_q95": float(sub["multiplier"].quantile(0.95)),
                    "delta_pct_q05": float(sub["delta_pct"].quantile(0.05)),
                    "delta_pct_q50": float(sub["delta_pct"].quantile(0.50)),
                    "delta_pct_q95": float(sub["delta_pct"].quantile(0.95)),
                }
            )
    out = pd.DataFrame(rows)
    val = validation_summary[validation_summary["method"].eq("overall_site_space_hgb")].copy()
    if not val.empty:
        val = val[["target_label", "target_key", "weighted_rmse", "rmse_reduction_vs_zero_pct", "rmse_reduction_vs_mean_pct"]]
        val["scope"] = "station_validation"
        out = pd.concat([out, val], ignore_index=True, sort=False)
    return out


def save_figure(summary: pd.DataFrame, validation_summary: pd.DataFrame, mesh_values: pd.DataFrame, station_values: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.9), constrained_layout=True)

    ax = axes[0, 0]
    val = validation_summary[validation_summary["method"].eq("overall_site_space_hgb")].copy()
    val["period_plot"] = val["target_key"].map({"sa03": 0.3, "sa10": 1.0, "sa30": 3.0})
    val = val.sort_values("period_plot")
    ax.bar(np.arange(len(val)), val["rmse_reduction_vs_zero_pct"], color="#4C78A8", label="vs zero")
    ax.scatter(np.arange(len(val)), val["rmse_reduction_vs_mean_pct"], color="#C44E52", zorder=3, label="vs mean")
    ax.axhline(0.0, color="0.25", lw=0.9)
    ax.set_xticks(np.arange(len(val)))
    ax.set_xticklabels(["0.3", "1.0", "3.0"])
    ax.set_xlabel("Anchor period (s)")
    ax.set_ylabel("RMSE reduction (%)")
    ax.set_title("A Spatial-block station validation")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    sub = summary[
        summary["scope"].eq("sampled_grid")
        & summary["probability_level"].eq("50y_10pct")
        & summary["period_s"].notna()
    ].copy()
    sub = sub.sort_values("period_s")
    ax.plot(sub["period_s"], sub["multiplier_q50"], color="#4C78A8", marker="o", lw=1.4, label="median")
    ax.fill_between(sub["period_s"], sub["multiplier_q05"], sub["multiplier_q95"], color="#4C78A8", alpha=0.20, label="5-95%")
    ax.axhline(1.0, color="0.25", lw=0.9)
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Multiplier")
    ax.set_title("B Sampled-grid correction spectrum")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    q10 = station_values[station_values["probability_level"].eq("50y_10pct")].copy()
    periods = sorted(q10["period_s"].unique())
    data = [q10[q10["period_s"].eq(period)]["delta_pct"].dropna().to_numpy() for period in periods]
    ax.boxplot(data, tick_labels=[f"{p:g}" for p in periods], showfliers=False)
    ax.axhline(0.0, color="0.25", lw=0.9)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Station delta (%)")
    ax.set_title("C Station-level official-UHS changes")

    ax = axes[1, 1]
    q10m = mesh_values[mesh_values["probability_level"].eq("50y_10pct")].copy()
    if "is_tail" in q10m:
        tail_mask = q10m["is_tail"].astype(bool)
        tail = q10m[tail_mask]
        nontail = q10m[~tail_mask]
    else:
        tail = q10m.iloc[0:0]
        nontail = q10m
    rows = []
    for label, frame in [("tail", tail), ("non-tail", nontail)]:
        for period_s, subp in frame.groupby("period_s"):
            rows.append({"label": label, "period_s": period_s, "median_delta_pct": subp["delta_pct"].median()})
    group = pd.DataFrame(rows)
    for label, color in [("tail", "#C44E52"), ("non-tail", "#4C78A8")]:
        g = group[group["label"].eq(label)].sort_values("period_s")
        if len(g):
            ax.plot(g["period_s"], g["median_delta_pct"], marker="o", lw=1.4, color=color, label=label)
    ax.axhline(0.0, color="0.25", lw=0.9)
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Median delta (%)")
    ax.set_title("D Residual-tail versus non-tail grid cells")
    ax.legend(frameon=False)

    fig.savefig(OUT_FIGURE, dpi=300)
    plt.close(fig)


def write_audit(
    stations: pd.DataFrame,
    mesh: pd.DataFrame,
    station_values: pd.DataFrame,
    mesh_values: pd.DataFrame,
    summary: pd.DataFrame,
    validation_summary: pd.DataFrame,
) -> None:
    val = validation_summary[validation_summary["method"].eq("overall_site_space_hgb")].copy()
    val_rows = []
    for _, row in val.iterrows():
        val_rows.append(
            f"- {row['target_label']}: weighted RMSE {row['weighted_rmse']:.4f} log10, "
            f"RMSE reduction {row['rmse_reduction_vs_zero_pct']:.1f}% versus zero correction and "
            f"{row['rmse_reduction_vs_mean_pct']:.1f}% versus national mean."
        )
    q10 = summary[
        summary["scope"].eq("sampled_grid")
        & summary["probability_level"].eq("50y_10pct")
        & summary["period_s"].notna()
    ].copy()
    lines = [
        "# All-period official-response-spectrum non-ergodic PSHA audit",
        "",
        "## Completed scope",
        "",
        f"- Station correction anchors: {len(stations):,} station-target rows across PGA, SA(0.3 s), SA(1.0 s), and SA(3.0 s).",
        f"- Response-spectrum anchor stations: {stations['siteid2'].nunique():,}.",
        f"- Sampled-grid cells: {len(mesh):,}.",
        f"- Official station response rows after correction: {len(station_values):,}.",
        f"- Official sampled-grid response rows after correction: {len(mesh_values):,}.",
        "- Official response-spectrum periods: 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, and 5.0 s.",
        "",
        "## Station validation",
        "",
        *val_rows,
        "",
        "## 50-year 10% sampled-grid multiplier spectrum",
        "",
        "| period_s | source | q05 | q50 | q95 | delta_pct_q50 |",
        "| ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for _, row in q10.sort_values("period_s").iterrows():
        lines.append(
            f"| {row['period_s']:.1f} | {row['correction_source']} | "
            f"{row['multiplier_q05']:.3f} | {row['multiplier_q50']:.3f} | {row['multiplier_q95']:.3f} | "
            f"{row['delta_pct_q50']:.1f} |"
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This audit propagates non-ergodic station terms through official J-SHIS response-spectrum map ordinates at the matched 250 m meshes. The 0.3, 1.0, and 3.0 s corrections are direct residual-model anchors. The 0.1 and 0.2 s corrections are clamped to the 0.3 s anchor; the 0.5 and 2.0 s corrections are log-period interpolations; the 5.0 s correction is clamped to the 3.0 s anchor. This closes an all-period official-map sensitivity calculation. It does not close a source-level official J-SHIS PSHA rerun because the complete source recurrence, period-dependent sigma, spatial correlation, and logic-tree implementation are not present in the local public-data package.",
            "",
            "## Output tables",
            "",
            f"- `{OUT_STATION_CORR.relative_to(ROOT)}`",
            f"- `{OUT_VALIDATION.relative_to(ROOT)}`",
            f"- `{OUT_VALIDATION_SUMMARY.relative_to(ROOT)}`",
            f"- `{OUT_MESH.relative_to(ROOT)}`",
            f"- `{OUT_STATION_VALUES.relative_to(ROOT)}`",
            f"- `{OUT_MESH_VALUES.relative_to(ROOT)}`",
            f"- `{OUT_SUMMARY.relative_to(ROOT)}`",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_to_supplement(paths: list[Path]) -> None:
    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    for path in paths:
        shutil.copy2(path, SUPPLEMENT_DIR / path.name)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    stations = load_station_corrections()
    validation, validation_summary = validate_surfaces(stations)
    mesh = build_mesh_surface(stations)

    station_anchor = build_station_anchor_wide(stations)
    mesh_anchor = build_mesh_anchor_wide(mesh)
    station_corrections = interpolate_corrections(station_anchor)
    mesh_corrections = interpolate_corrections(mesh_anchor)

    station_values = apply_response_maps(station_anchor, station_corrections, scope="station")
    mesh_values = apply_response_maps(mesh_anchor, mesh_corrections, scope="sampled_grid")
    for col in ["is_tail", "is_ne_hokkaido", "basin_period_proxy_1_3s", "nearest_station_distance_km", "nearest_siteid2", "nearest_site_code"]:
        if col in mesh_anchor.columns and col not in mesh_values.columns:
            mesh_values = mesh_values.merge(mesh_anchor[["meshcode250_str", col]], on="meshcode250_str", how="left", validate="many_to_one")
    summary = build_summary(station_values, mesh_values, validation_summary)

    stations.to_csv(OUT_STATION_CORR, index=False)
    validation.to_csv(OUT_VALIDATION, index=False)
    validation_summary.to_csv(OUT_VALIDATION_SUMMARY, index=False)
    mesh.to_csv(OUT_MESH, index=False)
    station_values.to_csv(OUT_STATION_VALUES, index=False)
    mesh_values.to_csv(OUT_MESH_VALUES, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    save_figure(summary, validation_summary, mesh_values, station_values)
    write_audit(stations, mesh, station_values, mesh_values, summary, validation_summary)
    copy_to_supplement(
        [
            Path(__file__),
            OUT_STATION_CORR,
            OUT_VALIDATION,
            OUT_VALIDATION_SUMMARY,
            OUT_MESH,
            OUT_STATION_VALUES,
            OUT_MESH_VALUES,
            OUT_SUMMARY,
            OUT_AUDIT,
        ]
    )

    q10 = summary[
        summary["scope"].eq("sampled_grid")
        & summary["probability_level"].eq("50y_10pct")
        & summary["period_s"].notna()
    ].sort_values("period_s")
    print(
        "Multiperiod official-map nonergodic PSHA audit complete: "
        f"{stations['siteid2'].nunique()} stations, {len(mesh)} sampled meshes, "
        f"{len(station_values)} station response rows, {len(mesh_values)} mesh response rows. "
        "50y10pct median multipliers by period: "
        + ", ".join(f"{row.period_s:g}s={row.multiplier_q50:.3f}" for row in q10.itertuples())
    )


if __name__ == "__main__":
    main()
