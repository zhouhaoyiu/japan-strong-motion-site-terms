#!/usr/bin/env python3
"""Transparent physical-proxy model for SA(3.0 s) station residuals."""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
ARTICLE = OUT / "cee_submission_latex_v0_8_english_article"
FIGURES = ARTICLE / "figures"
SUPPLEMENT = ARTICLE / "supplement"

PREDICTIONS = OUT / "jshis_station_model_ablation_predictions.csv"
SITE = OUT / "jshis_site_schema_v2024_sub1.csv"
BOREHOLE = SUPPLEMENT / "jshis_public_psha_10467_nearest_station_borehole_profiles.csv"

TARGET = "SA(3.0s) RotD50"
NUMERIC = [
    "log_d1400",
    "log_dbase",
    "log_avs30",
    "log_vs20",
    "log_borehole_vs20",
    "t_d1400_vs700_s",
    "basin_1_3s_window",
    "log_period_misfit_2s",
    "soft_upper20_thickness",
    "depth_to_vs_ge_500_m",
    "dist_vf_mf13_nejapan",
    "dist_vf_mf13_swjapan",
    "elevation",
    "sensor_depth_glminus",
    "lon",
    "lat",
]
CATEGORICAL = ["region"]


def region_name(lat: float, lon: float) -> str:
    if lat >= 41:
        return "Hokkaido"
    if lat >= 37:
        return "Tohoku"
    if lat >= 35 and lon >= 138:
        return "Kanto"
    if lat >= 34.5 and lon >= 135:
        return "Chubu_Kinki"
    if lat >= 33 and lon >= 130:
        return "Chugoku_Shikoku"
    return "Kyushu_Okinawa"


def weighted_rmse(y: np.ndarray, yhat: np.ndarray, w: np.ndarray) -> float:
    return float(np.sqrt(np.average((y - yhat) ** 2, weights=w)))


def make_encoder() -> OneHotEncoder:
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def make_model() -> Pipeline:
    return Pipeline(
        [
            (
                "prep",
                ColumnTransformer(
                    [
                        ("num", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), NUMERIC),
                        ("cat", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("onehot", make_encoder())]), CATEGORICAL),
                    ]
                ),
            ),
            ("ridge", RidgeCV(alphas=np.logspace(-3, 3, 25))),
        ]
    )


def load_data() -> pd.DataFrame:
    pred = pd.read_csv(PREDICTIONS)
    pred = pred[pred["split"].eq("spatial_block") & pred["target_label"].eq(TARGET)].copy()
    site = pd.read_csv(SITE).sort_values(["siteid2", "site_code"]).drop_duplicates("siteid2")
    borehole = pd.read_csv(BOREHOLE).sort_values("station_code").drop_duplicates("station_code")
    df = pred.merge(site, on=["siteid2", "site_code"], how="left").merge(
        borehole.add_prefix("borehole_"), left_on="site_code", right_on="borehole_station_code", how="left"
    )
    for col in [
        "lon",
        "lat",
        "d1400",
        "dbase",
        "avs30",
        "vs20",
        "dist_vf_mf13_nejapan",
        "dist_vf_mf13_swjapan",
        "elevation",
        "sensor_depth_glminus",
        "borehole_vs20_profile_mps",
        "borehole_low_vs_le_300_thickness_20m",
        "borehole_depth_to_vs_ge_500_m",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["region"] = [region_name(lat, lon) for lat, lon in zip(df["lat"], df["lon"])]
    df["log_d1400"] = np.log10(df["d1400"].clip(lower=1))
    df["log_dbase"] = np.log10(df["dbase"].clip(lower=1))
    df["log_avs30"] = np.log10(df["avs30"].clip(lower=1))
    df["log_vs20"] = np.log10(df["vs20"].clip(lower=1))
    df["log_borehole_vs20"] = np.log10(df["borehole_vs20_profile_mps"].clip(lower=1))
    df["t_d1400_vs700_s"] = 4.0 * df["d1400"] / 700.0
    df["basin_1_3s_window"] = ((df["t_d1400_vs700_s"] >= 1.0) & (df["t_d1400_vs700_s"] <= 3.0)).astype(float)
    df["log_period_misfit_2s"] = np.abs(np.log(df["t_d1400_vs700_s"].clip(lower=0.05) / 2.0))
    df["soft_upper20_thickness"] = df["borehole_low_vs_le_300_thickness_20m"]
    df["depth_to_vs_ge_500_m"] = df["borehole_depth_to_vs_ge_500_m"]
    return df


def fit_oof(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, Pipeline]:
    y = df["mean_residual"].astype(float).to_numpy()
    w = df["n_records"].astype(float).clip(lower=1).to_numpy()
    fold = df["fold"].astype(int).to_numpy()
    yhat = np.full(len(df), np.nan)
    for f in sorted(set(fold)):
        train = fold != f
        test = fold == f
        model = make_model()
        model.fit(df.loc[train, NUMERIC + CATEGORICAL], y[train], ridge__sample_weight=w[train])
        yhat[test] = model.predict(df.loc[test, NUMERIC + CATEGORICAL])
    zero = weighted_rmse(y, np.zeros_like(y), w)
    mean = np.average(y, weights=w)
    mean_rmse = weighted_rmse(y, np.full_like(y, mean), w)
    rmse = weighted_rmse(y, yhat, w)
    metrics = pd.DataFrame(
        [
            {
                "target_label": TARGET,
                "n_stations": len(df),
                "weighted_rmse_zero": zero,
                "weighted_rmse_national_mean": mean_rmse,
                "weighted_rmse_physical_proxy": rmse,
                "rmse_reduction_vs_zero_pct": (1.0 - rmse / zero) * 100.0,
                "rmse_reduction_vs_mean_pct": (1.0 - rmse / mean_rmse) * 100.0,
                "pearson_observed_predicted": float(np.corrcoef(y, yhat)[0, 1]),
            }
        ]
    )
    out = df[["siteid2", "site_code", "network_label", "n_records", "fold", "mean_residual"]].copy()
    out["physical_proxy_prediction"] = yhat
    out["physical_proxy_residual"] = out["mean_residual"] - yhat
    final_model = make_model()
    final_model.fit(df[NUMERIC + CATEGORICAL], y, ridge__sample_weight=w)
    return metrics, out, final_model


def coefficients(model: Pipeline) -> pd.DataFrame:
    prep = model.named_steps["prep"]
    num_names = NUMERIC
    cat_names = list(prep.named_transformers_["cat"].named_steps["onehot"].get_feature_names_out(CATEGORICAL))
    coefs = model.named_steps["ridge"].coef_
    return pd.DataFrame({"feature": num_names + cat_names, "standardized_coefficient": coefs}).sort_values(
        "standardized_coefficient", key=lambda s: s.abs(), ascending=False
    )


def plot(pred: pd.DataFrame, coef: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    ax = axes[0]
    ax.scatter(pred["mean_residual"], pred["physical_proxy_prediction"], s=9, alpha=0.35, edgecolors="none")
    lim = np.nanmax(np.abs(pred[["mean_residual", "physical_proxy_prediction"]].to_numpy())) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="#333333", linewidth=1.0)
    ax.axhline(0, color="#cccccc", linewidth=0.8)
    ax.axvline(0, color="#cccccc", linewidth=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Observed SA(3.0 s) station residual")
    ax.set_ylabel("Physical-proxy prediction")
    ax.set_title("A  Spatial-block held-out prediction")

    ax = axes[1]
    top = coef.head(12).iloc[::-1]
    colors = np.where(top["standardized_coefficient"] >= 0, "#C44E52", "#4C78A8")
    ax.barh(top["feature"], top["standardized_coefficient"], color=colors)
    ax.axvline(0, color="#333333", linewidth=0.9)
    ax.set_xlabel("Standardized ridge coefficient")
    ax.set_title("B  Physical and geographic proxy weights")
    fig.tight_layout()
    fig.savefig(FIGURES / "figure18_sa3_physical_proxy_model.png", dpi=300)
    fig.savefig(FIGURES / "figure18_sa3_physical_proxy_model.pdf")


def markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = df.round(4).astype(str).values.tolist()
    out = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    out.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(out)


def main() -> None:
    df = load_data()
    metrics, pred, model = fit_oof(df)
    coef = coefficients(model)
    plot(pred, coef)
    files = {
        "metrics": OUT / "jshis_sa3_physical_proxy_model_metrics.csv",
        "predictions": OUT / "jshis_sa3_physical_proxy_model_predictions.csv",
        "coefficients": OUT / "jshis_sa3_physical_proxy_model_coefficients.csv",
        "audit": OUT / "jshis_sa3_physical_proxy_model.md",
    }
    metrics.to_csv(files["metrics"], index=False)
    pred.to_csv(files["predictions"], index=False)
    coef.to_csv(files["coefficients"], index=False)
    files["audit"].write_text(
        "\n".join(
            [
                "# SA(3.0 s) physical-proxy station model",
                "",
                "This check fits a transparent ridge model to the MF2013 site-corrected SA(3.0 s) station residuals. Predictors are basin-depth proxies, shallow-velocity proxies, public borehole-profile proxies, volcanic-front distances, elevation, sensor depth, coordinates, and broad geographic region. No station identifier enters the model.",
                "",
                "## Spatial-block held-out performance",
                "",
                markdown_table(metrics),
                "",
                "## Largest standardized coefficients",
                "",
                markdown_table(coef.head(12)),
                "",
                "## Interpretation boundary",
                "",
                "The model is an independent physical and geographic proxy check. Calibrated three-dimensional basin-response simulation requires separate velocity and basin-geometry data. The proxy model shows that a simple interpretable model recovers part of the long-period station term outside the gradient-boosted station model.",
            ]
        )
        + "\n"
    )
    SUPPLEMENT.mkdir(parents=True, exist_ok=True)
    for path in files.values():
        shutil.copy2(path, SUPPLEMENT / path.name)
    print(files["audit"])


if __name__ == "__main__":
    main()
