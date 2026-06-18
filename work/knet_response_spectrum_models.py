#!/usr/bin/env python3
"""Leakage-safe models for K-NET response-spectrum targets."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_INPUT = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_record_response_spectrum_features.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)

TARGETS = [
    ("log10_psa_0.3s_5pct_horizontal", "SA(0.3s)"),
    ("log10_psa_1s_5pct_horizontal", "SA(1.0s)"),
    ("log10_psa_3s_5pct_horizontal", "SA(3.0s)"),
]


def load_frame(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"_id": str, "EventName": str})
    df = df[
        df["Magnitude"].notna()
        & df["Distance"].notna()
        & df["Depth_km"].notna()
        & df["StationCode"].notna()
        & df["EventName"].notna()
        & df["Latitude"].notna()
        & df["Longitude"].notna()
        & df["StationLat"].notna()
        & df["StationLong"].notna()
    ].copy()
    df["hyp_distance_km"] = np.sqrt(
        df["Distance"].astype(float) ** 2 + df["Depth_km"].astype(float) ** 2
    )
    df["log10_hyp_distance_km"] = np.log10(df["hyp_distance_km"] + 1.0)
    return df


def ridge_pipeline(numeric: list[str], categorical: list[str] | None = None) -> Pipeline:
    categorical = categorical or []
    transformers = [
        (
            "num",
            Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]),
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
    return Pipeline(
        [
            ("features", ColumnTransformer(transformers=transformers, remainder="drop")),
            ("model", Ridge(alpha=5.0)),
        ]
    )


def hgb_pipeline(numeric: list[str]) -> Pipeline:
    return Pipeline(
        [
            (
                "features",
                ColumnTransformer(
                    [("num", Pipeline([("impute", SimpleImputer(strategy="median"))]), numeric)],
                    remainder="drop",
                ),
            ),
            (
                "model",
                HistGradientBoostingRegressor(
                    max_iter=300,
                    learning_rate=0.05,
                    max_leaf_nodes=31,
                    l2_regularization=0.05,
                    random_state=20260605,
                ),
            ),
        ]
    )


def build_models() -> dict[str, Pipeline]:
    physical = ["Magnitude", "log10_hyp_distance_km", "Depth_km"]
    geo = physical + [
        "Latitude",
        "Longitude",
        "StationLat",
        "StationLong",
        "StationHeight_m",
    ]
    return {
        "ridge_physical": ridge_pipeline(physical),
        "ridge_site_terms": ridge_pipeline(geo, ["StationCode"]),
        "hgb_geo_path": hgb_pipeline(geo),
    }


def split_indices(df: pd.DataFrame, split_col: str) -> tuple[np.ndarray, np.ndarray]:
    train_mask = df[split_col].isin(["train", "dev"]).to_numpy()
    test_mask = (df[split_col] == "test").to_numpy()
    return np.where(train_mask)[0], np.where(test_mask)[0]


def evaluate(df: pd.DataFrame, split: str, target: str, label: str, model_name: str, model: Pipeline) -> dict:
    work = df[df[target].notna()].copy()
    train_idx, test_idx = split_indices(work, split)
    y = work[target].to_numpy()
    model.fit(work.iloc[train_idx], y[train_idx])
    pred = model.predict(work.iloc[test_idx])
    return {
        "split": split,
        "target": target,
        "target_label": label,
        "model": model_name,
        "train_rows": int(len(train_idx)),
        "test_rows": int(len(test_idx)),
        "mae": float(mean_absolute_error(y[test_idx], pred)),
        "r2": float(r2_score(y[test_idx], pred)),
    }


def add_reductions(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split, target), sub in results.groupby(["split", "target"]):
        base = sub[sub["model"] == "ridge_physical"].iloc[0]
        for _, row in sub.iterrows():
            item = row.to_dict()
            item["mae_reduction_vs_physical"] = base["mae"] - row["mae"]
            item["percent_mae_reduction_vs_physical"] = (
                (base["mae"] - row["mae"]) / base["mae"] * 100.0
            )
            rows.append(item)
    return pd.DataFrame(rows)


def plot_results(reductions: pd.DataFrame, output_dir: Path) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for split in ["event_split", "time_split", "station_split"]:
        sub = reductions[(reductions["split"] == split) & (reductions["model"] != "ridge_physical")]
        pivot = sub.pivot(index="target_label", columns="model", values="percent_mae_reduction_vs_physical")
        plt.figure(figsize=(6.7, 3.8))
        im = plt.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
        plt.yticks(np.arange(len(pivot.index)), pivot.index)
        plt.xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=20, ha="right")
        plt.colorbar(im, label="MAE reduction vs physical (%)")
        plt.title(f"Response-spectrum model gains: {split}")
        plt.tight_layout()
        plt.savefig(figures / f"knet_response_spectrum_reduction_heatmap_{split}.png", dpi=220)
        plt.close()


def write_report(reductions: pd.DataFrame, output_dir: Path) -> None:
    lines = ["# K-NET Response Spectrum Models\n\n"]
    lines.append("Targets are 5%-damped waveform-derived pseudo spectral accelerations.\n\n")
    for split in ["event_split", "time_split", "station_split"]:
        lines.append(f"## {split}\n")
        sub = reductions[reductions["split"] == split]
        for _, label in TARGETS:
            rows = sub[sub["target_label"] == label].sort_values("mae")
            best = rows.iloc[0]
            physical = rows[rows["model"] == "ridge_physical"].iloc[0]
            lines.append(
                f"- {label}: best={best['model']} MAE={best['mae']:.4f}, "
                f"R2={best['r2']:.4f}, gain vs physical={best['percent_mae_reduction_vs_physical']:.1f}% "
                f"(physical MAE={physical['mae']:.4f})\n"
            )
        lines.append("\n")
    lines.append("## Guardrails\n")
    lines.append("- These are PSA values from demeaned accelerograms; final publication should document baseline correction.\n")
    lines.append("- SA targets are stronger engineering-ground-motion endpoints than Fourier proxy bands.\n")
    (output_dir / "knet_response_spectrum_models.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = load_frame(args.input)
    models = build_models()
    rows = []
    for split in ["event_split", "time_split", "station_split"]:
        for target, label in TARGETS:
            for model_name, model in models.items():
                rows.append(evaluate(df, split, target, label, model_name, model))
    results = pd.DataFrame(rows)
    results.to_csv(args.output_dir / "knet_response_spectrum_model_metrics.csv", index=False)
    reductions = add_reductions(results)
    reductions.to_csv(args.output_dir / "knet_response_spectrum_model_reductions.csv", index=False)
    plot_results(reductions, args.output_dir)
    write_report(reductions, args.output_dir)


if __name__ == "__main__":
    main()
