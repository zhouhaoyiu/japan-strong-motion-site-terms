#!/usr/bin/env python3
"""Frequency-dependent station terms and relocation-safe site interpretation."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)
WAVEFORM_FILE = DEFAULT_OUTPUT_DIR / "knet_record_waveform_features.csv"
RESPONSE_FILE = DEFAULT_OUTPUT_DIR / "knet_record_response_spectrum_features.csv"

TARGETS = [
    ("log10_pga_horizontal_from_waveform_gal", "PGA H", "scalar"),
    ("log10_arias_proxy_horizontal", "Arias H", "energy"),
    ("log10_spectral_amp_0.2_0.5_hz_horizontal", "0.2-0.5 Hz", "fourier"),
    ("log10_spectral_amp_0.5_1_hz_horizontal", "0.5-1 Hz", "fourier"),
    ("log10_spectral_amp_1_2_hz_horizontal", "1-2 Hz", "fourier"),
    ("log10_spectral_amp_2_5_hz_horizontal", "2-5 Hz", "fourier"),
    ("log10_spectral_amp_5_10_hz_horizontal", "5-10 Hz", "fourier"),
    ("log10_spectral_amp_10_20_hz_horizontal", "10-20 Hz", "fourier"),
    ("log10_psa_0.3s_5pct_horizontal", "SA 0.3s", "response"),
    ("log10_psa_1s_5pct_horizontal", "SA 1.0s", "response"),
    ("log10_psa_3s_5pct_horizontal", "SA 3.0s", "response"),
]


def load_targets() -> pd.DataFrame:
    waveform = pd.read_csv(WAVEFORM_FILE, dtype={"_id": str, "EventName": str})
    response_cols = [
        "_id",
        "log10_psa_0.3s_5pct_horizontal",
        "log10_psa_1s_5pct_horizontal",
        "log10_psa_3s_5pct_horizontal",
    ]
    response = pd.read_csv(RESPONSE_FILE, usecols=response_cols, dtype={"_id": str})
    df = waveform.merge(response, on="_id", how="left")
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
    df["EventName"] = df["EventName"].astype(str)
    df["StationCode"] = df["StationCode"].astype(str)
    df["hyp_distance_km"] = np.sqrt(
        df["Distance"].astype(float) ** 2 + df["Depth_km"].astype(float) ** 2
    )
    df["log10_hyp_distance_km"] = np.log10(df["hyp_distance_km"] + 1.0)
    df["midpoint_lat"] = (df["Latitude"].astype(float) + df["StationLat"].astype(float)) / 2.0
    df["midpoint_lon"] = (df["Longitude"].astype(float) + df["StationLong"].astype(float)) / 2.0
    return df


def station_term_pipeline() -> Pipeline:
    numeric = [
        "Magnitude",
        "log10_hyp_distance_km",
        "Depth_km",
        "Latitude",
        "Longitude",
        "StationLat",
        "StationLong",
        "StationHeight_m",
        "midpoint_lat",
        "midpoint_lon",
    ]
    categorical = ["EventName", "StationCode"]
    features = ColumnTransformer(
        [
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
        sparse_threshold=1.0,
    )
    return Pipeline([("features", features), ("ridge", Ridge(alpha=3.0))])


def extract_station_terms(model: Pipeline, target: str, label: str, group: str, df: pd.DataFrame) -> pd.DataFrame:
    names = model.named_steps["features"].get_feature_names_out()
    coefs = model.named_steps["ridge"].coef_
    pattern = re.compile(r"^cat__StationCode_(.+)$")
    rows = []
    for name, coef in zip(names, coefs):
        match = pattern.match(name)
        if match:
            rows.append(
                {
                    "StationCode": match.group(1),
                    "target": target,
                    "target_label": label,
                    "target_group": group,
                    "regularized_station_term": float(coef),
                }
            )
    out = pd.DataFrame(rows)
    counts = df.groupby("StationCode").size().rename("target_n_records").reset_index()
    return out.merge(counts, on="StationCode", how="left")


def fit_all_station_terms(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    all_terms = []
    fit_rows = []
    for target, label, group in TARGETS:
        work = df[df[target].notna()].copy()
        model = station_term_pipeline()
        y = work[target].to_numpy(float)
        model.fit(work, y)
        pred = model.predict(work)
        residual = y - pred
        fit_rows.append(
            {
                "target": target,
                "target_label": label,
                "target_group": group,
                "rows": len(work),
                "residual_std": float(np.std(residual, ddof=1)),
                "residual_variance": float(np.var(residual, ddof=1)),
            }
        )
        all_terms.append(extract_station_terms(model, target, label, group, work))
        print(f"Fitted station terms for {label}", flush=True)
    terms = pd.concat(all_terms, ignore_index=True)
    terms.to_csv(output_dir / "knet_multitarget_station_terms.csv", index=False)
    pd.DataFrame(fit_rows).to_csv(
        output_dir / "knet_multitarget_station_term_fit_summary.csv", index=False
    )
    return terms


def merge_site_terms(output_dir: Path, min_records: int) -> pd.DataFrame:
    terms = pd.read_csv(output_dir / "knet_multitarget_station_terms.csv")
    site = pd.read_csv(
        output_dir / f"relocation_safe_knet_station_site_parameters_min{min_records}.csv"
    )
    site = site[site["safe_for_site_parameter_join"].fillna(False)].copy()
    for col in ["safe_vs10_proxy_mps", "safe_vs20_proxy_mps", "safe_vs30_proxy_mps"]:
        site[col] = pd.to_numeric(site[col], errors="coerce")
    site["log10_vs10"] = np.log10(site["safe_vs10_proxy_mps"])
    site["log10_vs30"] = np.log10(site["safe_vs30_proxy_mps"])
    site["log10_vs30_over_vs10"] = np.log10(
        site["safe_vs30_proxy_mps"] / site["safe_vs10_proxy_mps"]
    )
    merged = terms.merge(
        site[
            [
                "station_code",
                "soil_profile_source",
                "safe_vs10_proxy_mps",
                "safe_vs20_proxy_mps",
                "safe_vs30_proxy_mps",
                "log10_vs10",
                "log10_vs30",
                "log10_vs30_over_vs10",
                "safe_site_elevation_m",
                "regularized_term_log10_pga",
            ]
        ],
        left_on="StationCode",
        right_on="station_code",
        how="inner",
    )
    merged.to_csv(output_dir / "knet_multitarget_station_terms_with_site.csv", index=False)
    return merged


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, seed: int = 20260605, n_boot: int = 1000) -> tuple[float, float, float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    rho, pvalue = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boot = []
    n = len(x)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        val, _ = stats.spearmanr(x[idx], y[idx])
        if np.isfinite(val):
            boot.append(val)
    low, high = np.percentile(boot, [2.5, 97.5])
    return float(rho), float(pvalue), float(low), float(high)


def site_correlations(merged: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    features = [
        ("log10(VS10)", "log10_vs10"),
        ("log10(VS30)", "log10_vs30"),
        ("log10(VS30/VS10)", "log10_vs30_over_vs10"),
    ]
    rows = []
    target_order = {target: idx for idx, (target, _, _) in enumerate(TARGETS)}
    for (target, label, group), sub in merged.groupby(["target", "target_label", "target_group"]):
        for feature_label, feature_col in features:
            data = sub[[feature_col, "regularized_station_term"]].replace([np.inf, -np.inf], np.nan).dropna()
            rho, pvalue, low, high = bootstrap_spearman(
                data[feature_col].to_numpy(float),
                data["regularized_station_term"].to_numpy(float),
            )
            rows.append(
                {
                    "target": target,
                    "target_label": label,
                    "target_group": group,
                    "target_order": target_order[target],
                    "feature": feature_label,
                    "n": len(data),
                    "spearman_rho": rho,
                    "pvalue": pvalue,
                    "bootstrap_ci95_low": low,
                    "bootstrap_ci95_high": high,
                }
            )
    corr = pd.DataFrame(rows).sort_values(["target_order", "feature"])
    corr.to_csv(output_dir / "knet_multitarget_station_site_correlations.csv", index=False)
    return corr


def write_figures(merged: pd.DataFrame, corr: pd.DataFrame, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)

    heat = corr.pivot(index="target_label", columns="feature", values="spearman_rho")
    order_labels = [label for _, label, _ in TARGETS if label in heat.index]
    heat = heat.loc[order_labels]
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    im = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=25, ha="right")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Site-Profile Correlation with Frequency-Dependent Station Terms")
    fig.colorbar(im, ax=ax, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(figure_dir / "knet_multitarget_station_site_correlation_heatmap.png", dpi=180)
    plt.close(fig)

    selected = ["0.5-1 Hz", "1-2 Hz", "SA 1.0s"]
    fig, axes = plt.subplots(1, len(selected), figsize=(12.0, 4.0), sharey=True)
    for ax, label in zip(axes, selected):
        sub = merged[merged["target_label"] == label].dropna(
            subset=["safe_vs10_proxy_mps", "regularized_station_term"]
        )
        ax.scatter(
            sub["safe_vs10_proxy_mps"],
            sub["regularized_station_term"],
            s=np.clip(np.sqrt(sub["target_n_records"]) * 10, 20, 100),
            c="#3c6e9f",
            alpha=0.7,
            linewidth=0.2,
            edgecolor="white",
        )
        ax.axhline(0, color="#444444", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_title(label)
        ax.set_xlabel("VS10 proxy (m/s)")
    axes[0].set_ylabel("Regularized station term")
    fig.tight_layout()
    fig.savefig(figure_dir / "knet_frequency_station_terms_vs_vs10.png", dpi=180)
    plt.close(fig)


def write_report(merged: pd.DataFrame, corr: pd.DataFrame, output_dir: Path) -> None:
    vs10 = corr[corr["feature"] == "log10(VS10)"].sort_values("target_order")
    strongest = corr.reindex(corr["spearman_rho"].abs().sort_values(ascending=False).index).head(8)
    lines = ["# Frequency-Dependent Station Terms and Site Profiles\n\n"]
    lines.append("## Data used\n")
    lines.append(
        f"- Station-target rows after relocation-safe site merge: {len(merged):,}.\n"
    )
    lines.append(
        f"- Unique K-NET stations with safe site profiles: {merged['StationCode'].nunique():,}.\n"
    )
    lines.append(
        "- Station terms are regularized all-data diagnostic terms after source/path/geographic controls.\n\n"
    )
    lines.append("## VS10 correlation by target\n")
    for _, row in vs10.iterrows():
        lines.append(
            f"- {row['target_label']}: rho={row['spearman_rho']:.3f} "
            f"(95% CI {row['bootstrap_ci95_low']:.3f} to "
            f"{row['bootstrap_ci95_high']:.3f}).\n"
        )
    lines.append("\n## Strongest site-profile associations\n")
    for _, row in strongest.iterrows():
        lines.append(
            f"- {row['target_label']} vs {row['feature']}: rho={row['spearman_rho']:.3f} "
            f"(95% CI {row['bootstrap_ci95_low']:.3f} to {row['bootstrap_ci95_high']:.3f}).\n"
        )
    lines.append("\n## Interpretation\n")
    lines.append(
        "The site-profile signal is frequency dependent and remains modest. "
        "Shallow velocity proxies are more informative than VS30 alone for several "
        "targets, which is consistent with near-surface amplification being only one "
        "part of the nonergodic station structure.\n\n"
    )
    lines.append("## Outputs\n")
    lines.append("- `knet_multitarget_station_terms.csv`\n")
    lines.append("- `knet_multitarget_station_terms_with_site.csv`\n")
    lines.append("- `knet_multitarget_station_site_correlations.csv`\n")
    lines.append("- `figures/knet_multitarget_station_site_correlation_heatmap.png`\n")
    lines.append("- `figures/knet_frequency_station_terms_vs_vs10.png`\n\n")
    lines.append("## Guardrail\n")
    lines.append(
        "These station terms are diagnostic fixed-effect terms, not held-out predictions. "
        "They should support physical interpretation only alongside the leakage-safe "
        "model results and official/published GMM residual checks.\n"
    )
    (output_dir / "knet_multitarget_station_site_explanation.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-records", type=int, default=10)
    args = parser.parse_args()

    df = load_targets()
    terms_path = args.output_dir / "knet_multitarget_station_terms.csv"
    if not terms_path.exists():
        fit_all_station_terms(df, args.output_dir)
    merged = merge_site_terms(args.output_dir, args.min_records)
    corr = site_correlations(merged, args.output_dir)
    write_figures(merged, corr, args.output_dir)
    write_report(merged, corr, args.output_dir)


if __name__ == "__main__":
    main()
