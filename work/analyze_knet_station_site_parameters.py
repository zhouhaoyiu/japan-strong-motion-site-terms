#!/usr/bin/env python3
"""Relate relocation-safe K-NET station terms to site-profile proxies."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from sklearn.model_selection import KFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int = 2000,
    seed: int = 20260605,
) -> tuple[float, float, float, float]:
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


def prepare_table(output_dir: Path, min_records: int) -> pd.DataFrame:
    path = output_dir / f"relocation_safe_knet_station_site_parameters_min{min_records}.csv"
    frame = pd.read_csv(path)
    frame = frame[frame["safe_for_site_parameter_join"].fillna(False)].copy()
    for column in [
        "safe_vs10_proxy_mps",
        "safe_vs20_proxy_mps",
        "safe_vs30_proxy_mps",
        "safe_site_elevation_m",
        "regularized_term_log10_pga",
        "mean_log10_pga_residual",
        "adjusted_term_n_records",
        "median_distance_km",
        "median_magnitude",
        "local_latitude_median",
        "local_longitude_median",
    ]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["log10_vs10"] = np.log10(frame["safe_vs10_proxy_mps"])
    frame["log10_vs20"] = np.log10(frame["safe_vs20_proxy_mps"])
    frame["log10_vs30"] = np.log10(frame["safe_vs30_proxy_mps"])
    frame["log10_vs30_over_vs10"] = np.log10(
        frame["safe_vs30_proxy_mps"] / frame["safe_vs10_proxy_mps"]
    )
    frame["soft_vs30_lt_300"] = frame["safe_vs30_proxy_mps"].lt(300).astype(int)
    frame["very_soft_vs10_lt_200"] = frame["safe_vs10_proxy_mps"].lt(200).astype(int)
    frame["relocation_corrected_profile"] = frame["soil_profile_source"].str.startswith(
        "relocation_", na=False
    )
    frame.to_csv(output_dir / "knet_station_site_parameter_analysis_table.csv", index=False)
    return frame


def correlation_table(frame: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    features = [
        ("VS10 proxy", "safe_vs10_proxy_mps"),
        ("VS20 proxy", "safe_vs20_proxy_mps"),
        ("VS30 proxy", "safe_vs30_proxy_mps"),
        ("log10(VS10)", "log10_vs10"),
        ("log10(VS20)", "log10_vs20"),
        ("log10(VS30)", "log10_vs30"),
        ("log10(VS30/VS10)", "log10_vs30_over_vs10"),
        ("soft VS30<300", "soft_vs30_lt_300"),
        ("very soft VS10<200", "very_soft_vs10_lt_200"),
        ("site elevation", "safe_site_elevation_m"),
    ]
    targets = [
        ("regularized station term", "regularized_term_log10_pga"),
        ("mean PGA residual", "mean_log10_pga_residual"),
    ]
    rows = []
    for target_name, target_col in targets:
        for feature_name, feature_col in features:
            subset = frame[[target_col, feature_col]].dropna()
            if len(subset) < 20 or subset[feature_col].nunique() < 2:
                continue
            rho, pvalue, low, high = bootstrap_spearman(
                subset[feature_col].to_numpy(float),
                subset[target_col].to_numpy(float),
            )
            rows.append(
                {
                    "target": target_name,
                    "feature": feature_name,
                    "n": len(subset),
                    "spearman_rho": rho,
                    "pvalue": pvalue,
                    "bootstrap_ci95_low": low,
                    "bootstrap_ci95_high": high,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(output_dir / "knet_station_site_correlations.csv", index=False)
    return out


def fit_wls(
    frame: pd.DataFrame,
    target_col: str,
    feature_cols: list[str],
    label: str,
) -> dict[str, float | str | int]:
    cols = [target_col, *feature_cols, "adjusted_term_n_records"]
    data = frame[cols].replace([np.inf, -np.inf], np.nan).dropna()
    y = data[target_col].to_numpy(float)
    x = sm.add_constant(data[feature_cols], has_constant="add")
    weights = np.sqrt(data["adjusted_term_n_records"].to_numpy(float))
    model = sm.WLS(y, x, weights=weights).fit(cov_type="HC3")
    result: dict[str, float | str | int] = {
        "target": target_col,
        "model": label,
        "n": len(data),
        "weighted_r2": float(model.rsquared),
        "aic": float(model.aic),
        "features": ",".join(feature_cols),
    }
    for name in x.columns:
        result[f"coef_{name}"] = float(model.params[name])
        result[f"p_{name}"] = float(model.pvalues[name])
    return result


def cross_validated_r2(
    frame: pd.DataFrame,
    target_col: str,
    feature_cols: list[str],
    seed: int = 20260605,
) -> float:
    cols = [target_col, *feature_cols]
    data = frame[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 50:
        return float("nan")
    x = data[feature_cols].to_numpy(float)
    y = data[target_col].to_numpy(float)
    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    scores = cross_val_score(model, x, y, cv=cv, scoring="r2")
    return float(np.mean(scores))


def regression_table(frame: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    model_specs = [
        ("VS30 only", ["log10_vs30"]),
        ("VS10 and VS30", ["log10_vs10", "log10_vs30"]),
        ("site profile", ["log10_vs10", "log10_vs30", "log10_vs30_over_vs10"]),
        (
            "site plus regional controls",
            [
                "log10_vs10",
                "log10_vs30",
                "log10_vs30_over_vs10",
                "safe_site_elevation_m",
                "local_latitude_median",
                "local_longitude_median",
            ],
        ),
    ]
    targets = ["regularized_term_log10_pga", "mean_log10_pga_residual"]
    rows = []
    for target_col in targets:
        for label, cols in model_specs:
            row = fit_wls(frame, target_col, cols, label)
            row["cv_r2_unweighted_ridge"] = cross_validated_r2(frame, target_col, cols)
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(output_dir / "knet_station_site_regression_models.csv", index=False)
    return out


def binned_summary(frame: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    data = frame.dropna(subset=["safe_vs30_proxy_mps", "regularized_term_log10_pga"]).copy()
    data["vs30_bin"] = pd.cut(
        data["safe_vs30_proxy_mps"],
        bins=[0, 200, 300, 450, 760, np.inf],
        labels=["<200", "200-300", "300-450", "450-760", ">=760"],
    )
    summary = (
        data.groupby("vs30_bin", observed=True)
        .agg(
            n=("station_code", "count"),
            median_station_term=("regularized_term_log10_pga", "median"),
            mean_station_term=("regularized_term_log10_pga", "mean"),
            median_mean_pga_residual=("mean_log10_pga_residual", "median"),
            median_vs30=("safe_vs30_proxy_mps", "median"),
        )
        .reset_index()
    )
    summary.to_csv(output_dir / "knet_station_site_vs30_binned_summary.csv", index=False)
    return summary


def write_figures(frame: pd.DataFrame, binned: pd.DataFrame, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    data = frame.dropna(subset=["safe_vs30_proxy_mps", "regularized_term_log10_pga"]).copy()

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    colors = np.where(data["relocation_corrected_profile"], "#b45b5b", "#3c6e9f")
    sizes = np.clip(np.sqrt(data["adjusted_term_n_records"]) * 12, 25, 130)
    ax.scatter(
        data["safe_vs30_proxy_mps"],
        data["regularized_term_log10_pga"],
        s=sizes,
        c=colors,
        alpha=0.72,
        linewidth=0.2,
        edgecolor="white",
    )
    ax.axhline(0, color="#444444", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Relocation-safe VS30 proxy (m/s)")
    ax.set_ylabel("Regularized station term, log10 PGA")
    ax.set_title("Station Terms vs Relocation-Safe VS30")
    fig.tight_layout()
    fig.savefig(figure_dir / "knet_station_term_vs_safe_vs30.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.bar(binned["vs30_bin"].astype(str), binned["median_station_term"], color="#557f6f")
    ax.axhline(0, color="#444444", linewidth=0.8)
    for i, row in binned.iterrows():
        ax.text(i, row["median_station_term"], f"n={int(row['n'])}", ha="center", va="bottom" if row["median_station_term"] >= 0 else "top", fontsize=9)
    ax.set_xlabel("VS30 proxy bin (m/s)")
    ax.set_ylabel("Median regularized station term")
    ax.set_title("Binned Site-Profile Signal in Station Terms")
    fig.tight_layout()
    fig.savefig(figure_dir / "knet_station_term_vs30_bins.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    data2 = frame.dropna(subset=["safe_vs10_proxy_mps", "regularized_term_log10_pga"])
    ax.scatter(
        data2["safe_vs10_proxy_mps"],
        data2["regularized_term_log10_pga"],
        s=np.clip(np.sqrt(data2["adjusted_term_n_records"]) * 12, 25, 130),
        c="#7a6a9b",
        alpha=0.7,
        linewidth=0.2,
        edgecolor="white",
    )
    ax.axhline(0, color="#444444", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Relocation-safe VS10 proxy (m/s)")
    ax.set_ylabel("Regularized station term, log10 PGA")
    ax.set_title("Station Terms vs Shallow VS10")
    fig.tight_layout()
    fig.savefig(figure_dir / "knet_station_term_vs_safe_vs10.png", dpi=180)
    plt.close(fig)


def write_report(
    frame: pd.DataFrame,
    correlations: pd.DataFrame,
    regressions: pd.DataFrame,
    binned: pd.DataFrame,
    output_dir: Path,
    min_records: int,
) -> None:
    target_corr = correlations[
        correlations["target"].eq("regularized station term")
        & correlations["feature"].isin(["log10(VS10)", "log10(VS30)", "soft VS30<300", "very soft VS10<200"])
    ].copy()
    target_corr = target_corr.sort_values("spearman_rho")
    best_reg = regressions[
        regressions["target"].eq("regularized_term_log10_pga")
    ].sort_values("weighted_r2", ascending=False).head(1)
    lines = ["# K-NET Station Terms and Relocation-Safe Site Parameters\n\n"]
    lines.append("## Data used\n")
    lines.append(
        f"- K-NET stations with station-term support n_records >= {min_records}: {len(frame):,}.\n"
    )
    lines.append(
        f"- Relocation-corrected profiles: {int(frame['relocation_corrected_profile'].sum()):,}.\n"
    )
    lines.append(
        "- Target station term: regularized log10 PGA station term from source-path-site decomposition.\n\n"
    )
    lines.append("## Main statistical signal\n")
    for _, row in target_corr.iterrows():
        lines.append(
            f"- {row['feature']}: Spearman rho={row['spearman_rho']:.3f} "
            f"(95% bootstrap CI {row['bootstrap_ci95_low']:.3f} to "
            f"{row['bootstrap_ci95_high']:.3f}, n={int(row['n'])}).\n"
        )
    if not best_reg.empty:
        row = best_reg.iloc[0]
        lines.append(
            f"- Best weighted linear model among tested specifications: {row['model']}, "
            f"weighted R2={row['weighted_r2']:.3f}, 5-fold ridge CV R2="
            f"{row['cv_r2_unweighted_ridge']:.3f}.\n"
        )
    lines.append("\n## VS30 bins\n")
    for _, row in binned.iterrows():
        lines.append(
            f"- VS30 {row['vs30_bin']} m/s: n={int(row['n'])}, "
            f"median station term={row['median_station_term']:.3f} log10 units.\n"
        )
    lines.append("\n## Interpretation\n")
    lines.append(
        "Lower shallow S-wave velocities tend to correspond to more positive K-NET "
        "station terms, consistent with a site-amplification contribution. The "
        "relationship is real but not exhaustive: simple VS proxies explain only a "
        "limited fraction of station-term variance, so residual terms should still "
        "be interpreted as mixed site, path, basin, and metadata/model effects.\n\n"
    )
    lines.append("## Outputs\n")
    lines.append("- `knet_station_site_parameter_analysis_table.csv`\n")
    lines.append("- `knet_station_site_correlations.csv`\n")
    lines.append("- `knet_station_site_regression_models.csv`\n")
    lines.append("- `knet_station_site_vs30_binned_summary.csv`\n")
    lines.append("- `figures/knet_station_term_vs_safe_vs30.png`\n")
    lines.append("- `figures/knet_station_term_vs30_bins.png`\n")
    lines.append("- `figures/knet_station_term_vs_safe_vs10.png`\n\n")
    lines.append("## Guardrail\n")
    lines.append(
        "This is an interpretive station-level analysis, not proof that VS30 alone "
        "causes the station terms. It should be strengthened with official J-SHIS "
        "site_schema.tsv values and published Japan GMM residuals before submission.\n"
    )
    (output_dir / "knet_station_site_parameter_explanation.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-records", type=int, default=10)
    args = parser.parse_args()

    frame = prepare_table(args.output_dir, args.min_records)
    correlations = correlation_table(frame, args.output_dir)
    regressions = regression_table(frame, args.output_dir)
    binned = binned_summary(frame, args.output_dir)
    write_figures(frame, binned, args.output_dir)
    write_report(frame, correlations, regressions, binned, args.output_dir, args.min_records)


if __name__ == "__main__":
    main()
