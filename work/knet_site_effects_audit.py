#!/usr/bin/env python3
"""Initial K-NET audit for source-path-site strong-motion analysis.

This script is intentionally conservative:
- It reads only metadata from the BSON header collection.
- It counts accelerogram records without decoding all waveform arrays.
- It creates event-level, station-level, and time-based split columns.
- It fits a simple ergodic baseline for log10(PGA) and summarizes residual
  station/event terms for later non-ergodic modeling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from pathlib import Path
from typing import Iterable

import bson
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_KNET_DIR = Path("/Users/yojironoda/Downloads/s7rk7bj3zn-1/knet_1530")
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def iter_bson_documents(path: Path) -> Iterable[dict]:
    """Yield BSON documents from a MongoDB bson dump file."""
    with path.open("rb") as handle:
        while True:
            prefix = handle.read(4)
            if not prefix:
                break
            if len(prefix) != 4:
                raise ValueError(f"Truncated BSON length prefix in {path}")
            size = struct.unpack("<i", prefix)[0]
            if size < 5:
                raise ValueError(f"Invalid BSON document size {size} in {path}")
            payload = prefix + handle.read(size - 4)
            if len(payload) != size:
                raise ValueError(f"Truncated BSON document in {path}")
            yield bson.loads(payload)


def count_bson_documents(path: Path) -> tuple[int, list[int]]:
    """Count BSON documents quickly using length prefixes only."""
    count = 0
    first_sizes: list[int] = []
    with path.open("rb") as handle:
        while True:
            prefix = handle.read(4)
            if not prefix:
                break
            if len(prefix) != 4:
                raise ValueError(f"Truncated BSON length prefix in {path}")
            size = struct.unpack("<i", prefix)[0]
            if count < 5:
                first_sizes.append(size)
            handle.seek(size - 4, 1)
            count += 1
    return count, first_sizes


def stable_hash_split(value: str, train: int = 70, dev: int = 10) -> str:
    """Deterministic split from a string key."""
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()
    bucket = int(digest[:8], 16) % 100
    if bucket < train:
        return "train"
    if bucket < train + dev:
        return "dev"
    return "test"


def load_header_dataframe(knet_dir: Path) -> pd.DataFrame:
    rows = list(iter_bson_documents(knet_dir / "header.bson"))
    df = pd.DataFrame(rows)
    time_cols = ["OriginTime", "RecordTime", "pTime", "sTime"]
    for col in time_cols:
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], utc=True, errors="coerce")
    return df


def add_splits(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["event_split"] = out["EventName"].astype(str).map(stable_hash_split)
    out["station_split"] = out["StationCode"].astype(str).map(stable_hash_split)

    order = out[["EventName", "OriginTime"]].drop_duplicates().sort_values("OriginTime")
    n_events = len(order)
    train_cut = math.floor(n_events * 0.7)
    dev_cut = math.floor(n_events * 0.8)
    time_split_by_event = {}
    for idx, event in enumerate(order["EventName"].astype(str)):
        if idx < train_cut:
            split = "train"
        elif idx < dev_cut:
            split = "dev"
        else:
            split = "test"
        time_split_by_event[event] = split
    out["time_split"] = out["EventName"].astype(str).map(time_split_by_event)
    return out


def prepare_model_frame(df: pd.DataFrame) -> pd.DataFrame:
    keep = df.copy()
    keep = keep[
        keep["PGA_gal"].notna()
        & keep["Magnitude"].notna()
        & keep["Distance"].notna()
        & keep["Depth_km"].notna()
        & keep["StationCode"].notna()
        & keep["EventName"].notna()
    ].copy()
    keep = keep[keep["PGA_gal"] > 0]
    keep = keep[keep["Distance"] >= 0]
    keep["log10_pga_gal"] = np.log10(keep["PGA_gal"].astype(float))
    keep["log10_distance_km"] = np.log10(keep["Distance"].astype(float) + 1.0)
    keep["p_to_s_sec"] = (
        keep["sIndex"].astype(float) - keep["pIndex"].astype(float)
    ) / keep["SamplingRate_Hz"].astype(float)
    return keep


def make_ergodic_pipeline() -> Pipeline:
    numeric = ["Magnitude", "log10_distance_km", "Depth_km"]
    transformer = ColumnTransformer(
        transformers=[
            (
                "num",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="median")),
                        ("scale", StandardScaler()),
                    ]
                ),
                numeric,
            )
        ],
        remainder="drop",
    )
    return Pipeline([("features", transformer), ("ridge", Ridge(alpha=1.0))])


def make_nonergodic_pipeline() -> Pipeline:
    numeric = ["Magnitude", "log10_distance_km", "Depth_km"]
    categorical = ["StationCode", "EventName"]
    transformer = ColumnTransformer(
        transformers=[
            (
                "num",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="median")),
                        ("scale", StandardScaler()),
                    ]
                ),
                numeric,
            ),
            (
                "cat",
                OneHotEncoder(handle_unknown="ignore", min_frequency=2),
                categorical,
            ),
        ],
        remainder="drop",
    )
    return Pipeline([("features", transformer), ("ridge", Ridge(alpha=5.0))])


def evaluate_group_split(df: pd.DataFrame, group_col: str) -> dict:
    model_df = prepare_model_frame(df)
    x = model_df[["Magnitude", "log10_distance_km", "Depth_km"]]
    y = model_df["log10_pga_gal"].to_numpy()
    groups = model_df[group_col].astype(str).to_numpy()

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=20260605)
    train_idx, test_idx = next(splitter.split(x, y, groups=groups))
    model = make_ergodic_pipeline()
    model.fit(model_df.iloc[train_idx], y[train_idx])
    pred = model.predict(model_df.iloc[test_idx])
    return {
        "group": group_col,
        "train_rows": int(len(train_idx)),
        "test_rows": int(len(test_idx)),
        "test_mae_log10_pga": float(mean_absolute_error(y[test_idx], pred)),
        "test_r2_log10_pga": float(r2_score(y[test_idx], pred)),
    }


def fit_residual_terms(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    model_df = prepare_model_frame(df)
    y = model_df["log10_pga_gal"].to_numpy()

    ergodic = make_ergodic_pipeline()
    ergodic.fit(model_df, y)
    model_df["ergodic_pred_log10_pga"] = ergodic.predict(model_df)
    model_df["ergodic_residual_log10_pga"] = (
        model_df["log10_pga_gal"] - model_df["ergodic_pred_log10_pga"]
    )

    nonergodic = make_nonergodic_pipeline()
    nonergodic.fit(model_df, y)
    model_df["nonergodic_pred_log10_pga"] = nonergodic.predict(model_df)
    model_df["nonergodic_residual_log10_pga"] = (
        model_df["log10_pga_gal"] - model_df["nonergodic_pred_log10_pga"]
    )

    station_terms = (
        model_df.groupby("StationCode")
        .agg(
            n_records=("StationCode", "size"),
            station_latitude=("StationLat", "median"),
            station_longitude=("StationLong", "median"),
            station_height_m=("StationHeight_m", "median"),
            mean_log10_pga_residual=("ergodic_residual_log10_pga", "mean"),
            median_log10_pga_residual=("ergodic_residual_log10_pga", "median"),
            median_distance_km=("Distance", "median"),
            median_magnitude=("Magnitude", "median"),
        )
        .reset_index()
        .sort_values("mean_log10_pga_residual", ascending=False)
    )

    event_terms = (
        model_df.groupby("EventName")
        .agg(
            n_records=("EventName", "size"),
            origin_time=("OriginTime", "min"),
            magnitude=("Magnitude", "median"),
            depth_km=("Depth_km", "median"),
            latitude=("Latitude", "median"),
            longitude=("Longitude", "median"),
            mean_log10_pga_residual=("ergodic_residual_log10_pga", "mean"),
            median_log10_pga_residual=("ergodic_residual_log10_pga", "median"),
        )
        .reset_index()
        .sort_values("mean_log10_pga_residual", ascending=False)
    )

    metrics = {
        "ergodic_mae_log10_pga_all": float(
            mean_absolute_error(y, model_df["ergodic_pred_log10_pga"])
        ),
        "ergodic_r2_log10_pga_all": float(r2_score(y, model_df["ergodic_pred_log10_pga"])),
        "nonergodic_mae_log10_pga_all": float(
            mean_absolute_error(y, model_df["nonergodic_pred_log10_pga"])
        ),
        "nonergodic_r2_log10_pga_all": float(
            r2_score(y, model_df["nonergodic_pred_log10_pga"])
        ),
    }
    return station_terms, event_terms, metrics


def save_figures(df: pd.DataFrame, station_terms: pd.DataFrame, figures_dir: Path) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    model_df = prepare_model_frame(df)

    plt.figure(figsize=(7.2, 4.2))
    plt.hist(model_df["Magnitude"], bins=np.arange(2.5, 8.6, 0.25), color="#4c78a8")
    plt.xlabel("Magnitude")
    plt.ylabel("Record count")
    plt.title("K-NET record magnitude distribution")
    plt.tight_layout()
    plt.savefig(figures_dir / "knet_magnitude_histogram.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.8))
    scatter = plt.scatter(
        model_df["Distance"],
        model_df["PGA_gal"],
        c=model_df["Magnitude"],
        s=9,
        alpha=0.55,
        cmap="viridis",
        linewidths=0,
    )
    plt.yscale("log")
    plt.xlabel("Epicentral distance (km)")
    plt.ylabel("PGA (gal, log scale)")
    plt.title("PGA-distance scaling colored by magnitude")
    cbar = plt.colorbar(scatter)
    cbar.set_label("Magnitude")
    plt.tight_layout()
    plt.savefig(figures_dir / "knet_pga_distance_magnitude.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.5))
    coverage = model_df.groupby("StationCode").size().sort_values(ascending=False)
    plt.hist(coverage, bins=40, color="#59a14f")
    plt.xlabel("Records per station")
    plt.ylabel("Station count")
    plt.title("Station coverage")
    plt.tight_layout()
    plt.savefig(figures_dir / "knet_station_coverage.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 5.8))
    st = station_terms[station_terms["n_records"] >= 5].copy()
    vmax = np.nanpercentile(np.abs(st["mean_log10_pga_residual"]), 95)
    scatter = plt.scatter(
        st["station_longitude"],
        st["station_latitude"],
        c=st["mean_log10_pga_residual"],
        s=np.clip(st["n_records"], 8, 70),
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        edgecolors="k",
        linewidths=0.15,
        alpha=0.9,
    )
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title("Station mean residual after ergodic PGA baseline")
    cbar = plt.colorbar(scatter)
    cbar.set_label("Mean residual in log10(PGA)")
    plt.tight_layout()
    plt.savefig(figures_dir / "knet_station_residual_map.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.5))
    valid_ps = model_df["p_to_s_sec"].replace([np.inf, -np.inf], np.nan).dropna()
    plt.hist(valid_ps[(valid_ps > 0) & (valid_ps < 80)], bins=60, color="#f28e2b")
    plt.xlabel("S-P time (s)")
    plt.ylabel("Record count")
    plt.title("P/S separation from header picks")
    plt.tight_layout()
    plt.savefig(figures_dir / "knet_sp_time_histogram.png", dpi=220)
    plt.close()


def write_report(
    df: pd.DataFrame,
    station_terms: pd.DataFrame,
    event_terms: pd.DataFrame,
    metrics: dict,
    split_metrics: list[dict],
    accel_count: int,
    accel_sizes: list[int],
    output_dir: Path,
) -> None:
    model_df = prepare_model_frame(df)
    lines = []
    lines.append("# K-NET Initial Audit for Source-Path-Site Study\n")
    lines.append("## Data inventory\n")
    lines.append(f"- Header records: {len(df):,}\n")
    lines.append(f"- Accelerogram component records: {accel_count:,}\n")
    lines.append(f"- First accelerogram BSON document sizes: {accel_sizes}\n")
    lines.append(f"- Events: {df['EventName'].nunique():,}\n")
    lines.append(f"- Stations: {df['StationCode'].nunique():,}\n")
    lines.append(
        f"- Magnitude range: {model_df['Magnitude'].min():.2f} to {model_df['Magnitude'].max():.2f}\n"
    )
    lines.append(
        f"- Distance range: {model_df['Distance'].min():.3f} to {model_df['Distance'].max():.3f} km\n"
    )
    lines.append(
        f"- PGA range: {model_df['PGA_gal'].min():.3f} to {model_df['PGA_gal'].max():.3f} gal\n"
    )

    lines.append("\n## Split design\n")
    for col in ["event_split", "station_split", "time_split"]:
        counts = df[col].value_counts(dropna=False).to_dict()
        lines.append(f"- {col}: {json.dumps(counts, sort_keys=True)}\n")
    event_counts = df.drop_duplicates("EventName")["event_split"].value_counts().to_dict()
    time_event_counts = df.drop_duplicates("EventName")["time_split"].value_counts().to_dict()
    station_counts = df.drop_duplicates("StationCode")["station_split"].value_counts().to_dict()
    lines.append(
        f"- unique events by event_split: {json.dumps(event_counts, sort_keys=True)}\n"
    )
    lines.append(
        f"- unique events by time_split: {json.dumps(time_event_counts, sort_keys=True)}\n"
    )
    lines.append(
        f"- unique stations by station_split: {json.dumps(station_counts, sort_keys=True)}\n"
    )
    lines.append(
        "- Use event-level or time-based splits for predictive claims. Station terms fitted on all data are diagnostic, not out-of-sample performance claims.\n"
    )

    lines.append("\n## Baseline residual decomposition\n")
    for key, value in metrics.items():
        lines.append(f"- {key}: {value:.4f}\n")
    for item in split_metrics:
        lines.append(
            f"- Event/group holdout by {item['group']}: MAE={item['test_mae_log10_pga']:.4f}, "
            f"R2={item['test_r2_log10_pga']:.4f}, test rows={item['test_rows']:,}\n"
        )

    lines.append("\n## Largest positive station residuals after ergodic baseline\n")
    lines.append("- Listed stations require at least 10 records.\n")
    supported_station_terms = station_terms[station_terms["n_records"] >= 10]
    for _, row in supported_station_terms.head(10).iterrows():
        lines.append(
            f"- {row['StationCode']}: mean residual={row['mean_log10_pga_residual']:.3f}, "
            f"records={int(row['n_records'])}\n"
        )

    lines.append("\n## Largest positive event residuals after ergodic baseline\n")
    for _, row in event_terms.head(10).iterrows():
        lines.append(
            f"- {row['EventName']}: M={row['magnitude']:.2f}, depth={row['depth_km']:.1f} km, "
            f"mean residual={row['mean_log10_pga_residual']:.3f}, records={int(row['n_records'])}\n"
        )

    lines.append("\n## Next modeling steps\n")
    lines.append(
        "- Build a non-ergodic model with explicit source, path, and site terms and evaluate using event/time/station holdouts.\n"
    )
    lines.append(
        "- Add KiK-net surface/downhole pairs or the full NIED KiK-net archive to test whether station residuals correspond to physical site amplification.\n"
    )
    lines.append(
        "- Report predictive metrics only under leakage-safe splits; use all-data residual terms only as exploratory diagnostics.\n"
    )
    lines.append(
        "- Avoid claims about earthquake prediction. The defensible claim is strong-motion variability and site/path/source control.\n"
    )

    (output_dir / "knet_initial_audit.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--knet-dir", type=Path, default=DEFAULT_KNET_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = args.output_dir / "figures"

    header_df = load_header_dataframe(args.knet_dir)
    header_df = add_splits(header_df)
    header_df.to_csv(args.output_dir / "knet_header_index.csv", index=False)

    accel_count, accel_sizes = count_bson_documents(args.knet_dir / "accelerogram.bson")

    station_terms, event_terms, metrics = fit_residual_terms(header_df)
    station_terms.to_csv(args.output_dir / "knet_station_residual_terms.csv", index=False)
    event_terms.to_csv(args.output_dir / "knet_event_residual_terms.csv", index=False)

    split_metrics = [
        evaluate_group_split(header_df, "EventName"),
        evaluate_group_split(header_df, "StationCode"),
    ]
    pd.DataFrame(split_metrics).to_csv(
        args.output_dir / "knet_leakage_safe_baseline_metrics.csv", index=False
    )

    save_figures(header_df, station_terms, figures_dir)
    write_report(
        header_df,
        station_terms,
        event_terms,
        metrics,
        split_metrics,
        accel_count,
        accel_sizes,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
