#!/usr/bin/env python3
"""Regularized source-path-site decomposition for K-NET log10(PGA).

This is a diagnostic decomposition, not the final publication model. It uses
regularized fixed effects to quantify how much event and station repeatability
remain after a physical magnitude-distance-depth baseline.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Geod
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_header_index.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def load_frame(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["OriginTime", "RecordTime", "pTime", "sTime"])
    df = df[
        df["PGA_gal"].notna()
        & df["Magnitude"].notna()
        & df["Distance"].notna()
        & df["Depth_km"].notna()
        & df["StationCode"].notna()
        & df["EventName"].notna()
        & df["Latitude"].notna()
        & df["Longitude"].notna()
        & df["StationLat"].notna()
        & df["StationLong"].notna()
        & (df["PGA_gal"] > 0)
        & (df["Distance"] >= 0)
    ].copy()
    df["EventName"] = df["EventName"].astype(str)
    df["StationCode"] = df["StationCode"].astype(str)
    df["log10_pga_gal"] = np.log10(df["PGA_gal"].astype(float))
    df["hyp_distance_km"] = np.sqrt(
        df["Distance"].astype(float) ** 2 + df["Depth_km"].astype(float) ** 2
    )
    df["log10_distance_km"] = np.log10(df["Distance"].astype(float) + 1.0)
    df["log10_hyp_distance_km"] = np.log10(df["hyp_distance_km"] + 1.0)

    geod = Geod(ellps="WGS84")
    az12, az21, dist_m = geod.inv(
        df["Longitude"].to_numpy(),
        df["Latitude"].to_numpy(),
        df["StationLong"].to_numpy(),
        df["StationLat"].to_numpy(),
    )
    df["geodesic_distance_km"] = dist_m / 1000.0
    df["path_azimuth_deg"] = az12
    radians = np.deg2rad(az12)
    df["path_azimuth_sin"] = np.sin(radians)
    df["path_azimuth_cos"] = np.cos(radians)
    df["midpoint_lat"] = (df["Latitude"].astype(float) + df["StationLat"].astype(float)) / 2.0
    df["midpoint_lon"] = (df["Longitude"].astype(float) + df["StationLong"].astype(float)) / 2.0
    return df


def build_pipeline(numeric: list[str], categorical: list[str]) -> Pipeline:
    transformers = [
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
    ]
    if categorical:
        transformers.append(
            (
                "cat",
                OneHotEncoder(handle_unknown="ignore", min_frequency=2),
                categorical,
            )
        )
    features = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=1.0)
    return Pipeline([("features", features), ("ridge", Ridge(alpha=3.0))])


def evaluate_fit(name: str, model: Pipeline, df: pd.DataFrame) -> tuple[dict, np.ndarray]:
    y = df["log10_pga_gal"].to_numpy()
    model.fit(df, y)
    pred = model.predict(df)
    residual = y - pred
    return (
        {
            "model": name,
            "rows": len(df),
            "mae_log10_pga": float(mean_absolute_error(y, pred)),
            "r2_log10_pga": float(r2_score(y, pred)),
            "residual_variance": float(np.var(residual, ddof=1)),
            "residual_std": float(np.std(residual, ddof=1)),
        },
        residual,
    )


def feature_names(model: Pipeline) -> np.ndarray:
    transformer = model.named_steps["features"]
    return transformer.get_feature_names_out()


def extract_categorical_terms(
    model: Pipeline, prefix: str, output_col: str
) -> pd.DataFrame:
    names = feature_names(model)
    coefs = model.named_steps["ridge"].coef_
    rows = []
    pattern = re.compile(rf"^cat__{re.escape(prefix)}_(.+)$")
    for name, coef in zip(names, coefs):
        match = pattern.match(name)
        if match:
            rows.append({output_col: match.group(1), "regularized_term_log10_pga": float(coef)})
    return pd.DataFrame(rows)


def save_figures(summary: pd.DataFrame, station_terms: pd.DataFrame, df: pd.DataFrame, output_dir: Path) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    baseline_var = summary.loc[summary["model"] == "physical", "residual_variance"].iloc[0]
    plot_df = summary.copy()
    plot_df["variance_reduction_vs_physical_pct"] = (
        (baseline_var - plot_df["residual_variance"]) / baseline_var * 100.0
    )
    plt.figure(figsize=(8.4, 4.8))
    plt.bar(plot_df["model"], plot_df["variance_reduction_vs_physical_pct"], color="#4c78a8")
    plt.ylabel("Residual variance reduction vs physical baseline (%)")
    plt.title("Regularized source-path-site decomposition")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(figures / "knet_source_path_site_variance_reduction.png", dpi=220)
    plt.close()

    if len(station_terms):
        station_geo = (
            df.groupby("StationCode")
            .agg(
                station_latitude=("StationLat", "median"),
                station_longitude=("StationLong", "median"),
            )
            .reset_index()
        )
        st = station_terms.merge(station_geo, on="StationCode", how="left")
        st = st[st["n_records"] >= 10]
        vmax = np.nanpercentile(np.abs(st["regularized_term_log10_pga"]), 95)
        plt.figure(figsize=(7.2, 5.8))
        scatter = plt.scatter(
            st["station_longitude"],
            st["station_latitude"],
            c=st["regularized_term_log10_pga"],
            s=np.clip(st["n_records"], 8, 70),
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            edgecolors="k",
            linewidths=0.15,
        )
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title("Adjusted station terms after source-path controls")
        cbar = plt.colorbar(scatter)
        cbar.set_label("Regularized station term in log10(PGA)")
        plt.tight_layout()
        plt.savefig(figures / "knet_adjusted_station_terms_map.png", dpi=220)
        plt.close()


def write_report(summary: pd.DataFrame, station_terms: pd.DataFrame, event_terms: pd.DataFrame, output_dir: Path) -> None:
    physical_var = summary.loc[summary["model"] == "physical", "residual_variance"].iloc[0]
    lines = ["# K-NET Regularized Source-Path-Site Decomposition\n\n"]
    lines.append("## Model ladder\n")
    for _, row in summary.iterrows():
        reduction = (physical_var - row["residual_variance"]) / physical_var * 100.0
        lines.append(
            f"- {row['model']}: MAE={row['mae_log10_pga']:.4f}, R2={row['r2_log10_pga']:.4f}, "
            f"residual std={row['residual_std']:.4f}, variance reduction vs physical={reduction:.1f}%\n"
        )

    lines.append("\n## Strong positive adjusted station terms\n")
    for _, row in station_terms.sort_values("regularized_term_log10_pga", ascending=False).head(12).iterrows():
        lines.append(f"- {row['StationCode']}: station term={row['regularized_term_log10_pga']:.3f}\n")

    lines.append("\n## Strong positive adjusted event terms\n")
    for _, row in event_terms.sort_values("regularized_term_log10_pga", ascending=False).head(12).iterrows():
        lines.append(f"- {row['EventName']}: event term={row['regularized_term_log10_pga']:.3f}\n")

    lines.append("\n## Interpretation guardrails\n")
    lines.append(
        "- This is an all-data diagnostic decomposition using regularized fixed effects. It identifies repeatable event and station structure but is not an out-of-sample predictive claim.\n"
    )
    lines.append(
        "- Predictive claims must use event/time/station holdouts from the separate leakage-safe model comparison.\n"
    )
    lines.append(
        "- Publication-grade source-path-site interpretation requires official station metadata and multi-event KiK-net surface/downhole validation.\n"
    )
    (output_dir / "knet_source_path_site_decomposition.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_frame(args.input)
    physical = ["Magnitude", "log10_hyp_distance_km", "Depth_km"]
    path_geo = physical + [
        "Latitude",
        "Longitude",
        "StationLat",
        "StationLong",
        "StationHeight_m",
        "path_azimuth_sin",
        "path_azimuth_cos",
        "midpoint_lat",
        "midpoint_lon",
    ]
    models = {
        "physical": build_pipeline(physical, []),
        "source_event_terms": build_pipeline(physical, ["EventName"]),
        "site_station_terms": build_pipeline(path_geo, ["StationCode"]),
        "source_site_terms": build_pipeline(path_geo, ["EventName", "StationCode"]),
    }

    rows = []
    residuals = {}
    for name, model in models.items():
        row, residual = evaluate_fit(name, model, df)
        rows.append(row)
        residuals[name] = residual

    summary = pd.DataFrame(rows)
    summary.to_csv(args.output_dir / "knet_source_path_site_variance_decomposition.csv", index=False)

    full_model = models["source_site_terms"]
    station_terms = extract_categorical_terms(full_model, "StationCode", "StationCode")
    event_terms = extract_categorical_terms(full_model, "EventName", "EventName")
    station_counts = df.groupby("StationCode").size().rename("n_records").reset_index()
    event_counts = df.groupby("EventName").size().rename("n_records").reset_index()
    station_terms = station_terms.merge(station_counts, on="StationCode", how="left")
    event_terms = event_terms.merge(event_counts, on="EventName", how="left")
    station_terms.to_csv(args.output_dir / "knet_adjusted_station_terms.csv", index=False)
    event_terms.to_csv(args.output_dir / "knet_adjusted_event_terms.csv", index=False)

    save_figures(summary, station_terms, df, args.output_dir)
    write_report(summary, station_terms, event_terms, args.output_dir)


if __name__ == "__main__":
    main()
