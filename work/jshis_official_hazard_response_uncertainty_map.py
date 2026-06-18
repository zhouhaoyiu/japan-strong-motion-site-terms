#!/usr/bin/env python3
"""Build uncertainty tables and spatial figure for the official J-SHIS check."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OFFICIAL_VALUES = PROJECT_ROOT / "outputs" / "jshis_official_response_sa3_station_values.csv"
DEFAULT_IMPACT = PROJECT_ROOT / "outputs" / "jshis_station_hazard_impact_by_station.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures"

BOOTSTRAP_SEED = 20260612
BOOTSTRAP_REPS = 2000
MAIN_LEVEL = "50y_10pct"


LABELS = {
    "en": {
        "title": "Spatial stability and uncertainty of official SA(3.0 s) response-spectrum check",
        "panel_a": "A Spatial pattern of station multipliers",
        "panel_b": "B Spatial-block fold medians",
        "panel_c": "C Distribution shift at 50-year 10%",
        "panel_d": "D Bootstrap intervals for distribution quantiles",
        "multiplier": "SA(3.0 s) station multiplier",
        "longitude": "Longitude",
        "latitude": "Latitude",
        "fold": "Spatial-block fold",
        "sa3": "SA(3.0 s), 50-year 10% (g)",
        "official": "Official",
        "corrected": "Station-corrected",
        "stations": "Number of stations",
        "quantile": "Distribution quantile",
    },
    "zh": {
        "title": "官方 SA(3.0 s) 响应谱核验的空间稳定性和不确定性",
        "panel_a": "A 台站倍率空间分布",
        "panel_b": "B 空间分块 fold 中位数",
        "panel_c": "C 50 年 10% 谱值分布变化",
        "panel_d": "D 分布分位 bootstrap 区间",
        "multiplier": "SA(3.0 s) 台站谱值倍率",
        "longitude": "经度",
        "latitude": "纬度",
        "fold": "空间分块 fold",
        "sa3": "SA(3.0 s)，50 年 10% (g)",
        "official": "官方",
        "corrected": "台站校正",
        "stations": "台站数",
        "quantile": "分布分位",
    },
}


def load_frame(official_values: Path, impact_path: Path) -> pd.DataFrame:
    official = pd.read_csv(official_values)
    impact = pd.read_csv(impact_path, usecols=["siteid2", "site_code", "fold", "sa3_uhs_multiplier"])
    frame = official.merge(impact, on=["siteid2", "site_code"], how="left", validate="many_to_one")
    if "sa3_uhs_multiplier" in frame.columns:
        frame["station_multiplier"] = frame["station_multiplier"].fillna(frame["sa3_uhs_multiplier"])
    frame = frame.drop(columns=[col for col in ["sa3_uhs_multiplier"] if col in frame.columns])
    return frame


def summarize_folds(frame: pd.DataFrame) -> pd.DataFrame:
    q10 = frame[frame["probability_level"].eq(MAIN_LEVEL)].copy()
    rows = []
    for fold, sub in q10.groupby("fold", sort=True):
        rows.append(
            {
                "fold": int(fold),
                "n_stations": len(sub),
                "station_multiplier_q05": sub["station_multiplier"].quantile(0.05),
                "station_multiplier_q50": sub["station_multiplier"].quantile(0.50),
                "station_multiplier_q95": sub["station_multiplier"].quantile(0.95),
                "official_sa3_g_q50": sub["official_sa3_g"].quantile(0.50),
                "corrected_sa3_g_q50": sub["corrected_sa3_g"].quantile(0.50),
                "official_sa3_g_q05": sub["official_sa3_g"].quantile(0.05),
                "official_sa3_g_q95": sub["official_sa3_g"].quantile(0.95),
                "corrected_sa3_g_q05": sub["corrected_sa3_g"].quantile(0.05),
                "corrected_sa3_g_q95": sub["corrected_sa3_g"].quantile(0.95),
                "delta_pct_q50": sub["delta_pct"].quantile(0.50),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_quantile_ci(frame: pd.DataFrame, reps: int, seed: int) -> pd.DataFrame:
    q10 = frame[frame["probability_level"].eq(MAIN_LEVEL)].copy().reset_index(drop=True)
    rng = np.random.default_rng(seed)
    metrics = {
        "station_multiplier": q10["station_multiplier"].to_numpy(dtype=float),
        "official_sa3_g": q10["official_sa3_g"].to_numpy(dtype=float),
        "corrected_sa3_g": q10["corrected_sa3_g"].to_numpy(dtype=float),
        "delta_pct": q10["delta_pct"].to_numpy(dtype=float),
    }
    quantiles = [0.05, 0.50, 0.95]
    rows = []
    n = len(q10)
    for metric, values in metrics.items():
        boot = np.empty((reps, len(quantiles)), dtype=float)
        for idx in range(reps):
            sample = values[rng.integers(0, n, size=n)]
            boot[idx, :] = np.quantile(sample, quantiles)
        observed = np.quantile(values, quantiles)
        low = np.quantile(boot, 0.025, axis=0)
        high = np.quantile(boot, 0.975, axis=0)
        for q, obs, lo, hi in zip(quantiles, observed, low, high, strict=True):
            rows.append(
                {
                    "probability_level": MAIN_LEVEL,
                    "metric": metric,
                    "distribution_quantile": q,
                    "observed": obs,
                    "bootstrap_ci_low": lo,
                    "bootstrap_ci_high": hi,
                    "n_stations": n,
                    "bootstrap_reps": reps,
                    "seed": seed,
                }
            )
    return pd.DataFrame(rows)


def write_report(frame: pd.DataFrame, fold_summary: pd.DataFrame, boot: pd.DataFrame, output_path: Path) -> None:
    def markdown_table(data: pd.DataFrame) -> str:
        columns = list(data.columns)
        rows = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
        for _, item in data.iterrows():
            vals = []
            for col in columns:
                value = item[col]
                if isinstance(value, (float, np.floating)):
                    vals.append(f"{float(value):.6g}")
                else:
                    vals.append(str(value))
            rows.append("| " + " | ".join(vals) + " |")
        return "\n".join(rows)

    q10 = frame[frame["probability_level"].eq(MAIN_LEVEL)].copy()
    mult = q10["station_multiplier"]
    official = q10["official_sa3_g"]
    corrected = q10["corrected_sa3_g"]
    fold_min = fold_summary["corrected_sa3_g_q50"].min()
    fold_max = fold_summary["corrected_sa3_g_q50"].max()
    lines = [
        "# Official J-SHIS Response-Spectrum Uncertainty and Spatial Check",
        "",
        f"- Matched stations at {MAIN_LEVEL}: {len(q10)}",
        f"- Station multiplier q05/q50/q95: {mult.quantile(0.05):.3f} / {mult.quantile(0.50):.3f} / {mult.quantile(0.95):.3f}",
        f"- Official SA(3.0 s) q05/q50/q95: {official.quantile(0.05):.3f} / {official.quantile(0.50):.3f} / {official.quantile(0.95):.3f} g",
        f"- Corrected SA(3.0 s) q05/q50/q95: {corrected.quantile(0.05):.3f} / {corrected.quantile(0.50):.3f} / {corrected.quantile(0.95):.3f} g",
        f"- Spatial-block corrected median range across folds: {fold_min:.3f}--{fold_max:.3f} g",
        "",
        "This check uses station-level multipliers from the spatial-block held-out residual model. It adjusts only the SA(3.0 s) ordinate from the official J-SHIS response-spectrum map. It is not a full PSHA rerun.",
        "",
        "## Bootstrap confidence intervals",
        "",
        markdown_table(boot),
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def setup_style(language: str) -> None:
    fonts = ["DejaVu Sans"]
    if language == "zh":
        fonts = ["Hiragino Sans GB", "Songti SC", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": fonts,
            "axes.unicode_minus": False,
            "font.size": 9.2,
            "axes.titlesize": 10.8,
            "axes.labelsize": 9.3,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "legend.fontsize": 7.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def try_add_map_features(ax: plt.Axes) -> bool:
    try:
        import cartopy.crs as ccrs  # noqa: F401

        ax.coastlines(resolution="10m", linewidth=0.45, color="#555555")
        gl = ax.gridlines(draw_labels=True, linewidth=0.25, color="#D9D9D9", alpha=0.8)
        gl.top_labels = False
        gl.right_labels = False
        return True
    except Exception:
        return False


def build_figure(frame: pd.DataFrame, fold_summary: pd.DataFrame, boot: pd.DataFrame, figure_dir: Path, stem: str, language: str) -> None:
    labels = LABELS[language]
    setup_style(language)
    q10 = frame[frame["probability_level"].eq(MAIN_LEVEL)].copy()
    figure_dir.mkdir(parents=True, exist_ok=True)

    try:
        import cartopy.crs as ccrs

        projection = ccrs.PlateCarree()
        fig = plt.figure(figsize=(12.3, 8.7))
        grid = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.30)
        ax_a = fig.add_subplot(grid[0, 0], projection=projection)
        ax_a.set_extent([128.0, 146.5, 30.0, 46.5], crs=projection)
        try_add_map_features(ax_a)
        scatter_transform = projection
    except Exception:
        fig = plt.figure(figsize=(12.3, 8.7))
        grid = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.30)
        ax_a = fig.add_subplot(grid[0, 0])
        ax_a.set_xlim(128.0, 146.5)
        ax_a.set_ylim(30.0, 46.5)
        ax_a.grid(color="#E6E6E6")
        scatter_transform = None

    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    norm = TwoSlopeNorm(vmin=0.20, vcenter=1.0, vmax=1.55)
    scatter_kwargs = {
        "x": q10["lon"],
        "y": q10["lat"],
        "c": q10["station_multiplier"],
        "s": 9,
        "cmap": "coolwarm",
        "norm": norm,
        "alpha": 0.86,
        "linewidths": 0,
    }
    if scatter_transform is not None:
        scatter_kwargs["transform"] = scatter_transform
    sc = ax_a.scatter(**scatter_kwargs)
    ax_a.set_title(labels["panel_a"])
    ax_a.set_xlabel(labels["longitude"])
    ax_a.set_ylabel(labels["latitude"])
    cbar = fig.colorbar(sc, ax=ax_a, shrink=0.82, pad=0.03)
    cbar.set_label(labels["multiplier"])

    x = fold_summary["fold"].to_numpy(dtype=int)
    width = 0.32
    ax_b.bar(x - width / 2, fold_summary["official_sa3_g_q50"], width=width, color="#777777", label=labels["official"])
    ax_b.bar(x + width / 2, fold_summary["corrected_sa3_g_q50"], width=width, color="#C44E52", label=labels["corrected"])
    for _, row in fold_summary.iterrows():
        fold = int(row["fold"])
        ax_b.plot([fold - width / 2, fold - width / 2], [row["official_sa3_g_q05"], row["official_sa3_g_q95"]], color="#444444", lw=1.0)
        ax_b.plot([fold + width / 2, fold + width / 2], [row["corrected_sa3_g_q05"], row["corrected_sa3_g_q95"]], color="#8B1A1A", lw=1.0)
    ax_b.set_title(labels["panel_b"])
    ax_b.set_xlabel(labels["fold"])
    ax_b.set_ylabel(labels["sa3"])
    ax_b.set_xticks(x)
    ax_b.grid(axis="y", color="#E6E6E6")
    ax_b.legend(frameon=False)

    bins = np.linspace(0.0, 0.35, 36)
    ax_c.hist(q10["official_sa3_g"], bins=bins, color="#777777", alpha=0.45, label=labels["official"])
    ax_c.hist(q10["corrected_sa3_g"], bins=bins, color="#C44E52", alpha=0.55, label=labels["corrected"])
    ax_c.set_title(labels["panel_c"])
    ax_c.set_xlabel(labels["sa3"])
    ax_c.set_ylabel(labels["stations"])
    ax_c.grid(axis="y", color="#E6E6E6")
    ax_c.legend(frameon=False)

    quantile_labels = ["5%", "50%", "95%"]
    xpos = np.arange(len(quantile_labels))
    for metric, offset, color, marker, name in [
        ("official_sa3_g", -0.08, "#555555", "o", labels["official"]),
        ("corrected_sa3_g", 0.08, "#C44E52", "s", labels["corrected"]),
    ]:
        sub = boot[boot["metric"].eq(metric)].sort_values("distribution_quantile")
        observed = sub["observed"].to_numpy(dtype=float)
        low = sub["bootstrap_ci_low"].to_numpy(dtype=float)
        high = sub["bootstrap_ci_high"].to_numpy(dtype=float)
        yerr = np.vstack([observed - low, high - observed])
        ax_d.errorbar(xpos + offset, observed, yerr=yerr, color=color, marker=marker, lw=1.2, capsize=3, label=name)
    ax_d.set_xticks(xpos)
    ax_d.set_xticklabels(quantile_labels)
    ax_d.set_title(labels["panel_d"])
    ax_d.set_xlabel(labels["quantile"])
    ax_d.set_ylabel(labels["sa3"])
    ax_d.grid(axis="y", color="#E6E6E6")
    ax_d.legend(frameon=False)

    fig.suptitle(labels["title"], fontsize=13.0)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official-values", type=Path, default=DEFAULT_OFFICIAL_VALUES)
    parser.add_argument("--impact", type=Path, default=DEFAULT_IMPACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--figure-stem", default="cee_fig9_spatial_uncertainty")
    parser.add_argument("--language", choices=["en", "zh"], default="en")
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = load_frame(args.official_values, args.impact)
    fold_summary = summarize_folds(frame)
    boot = bootstrap_quantile_ci(frame, args.bootstrap_reps, BOOTSTRAP_SEED)

    fold_summary.to_csv(args.output_dir / "jshis_official_hazard_response_fold_summary.csv", index=False)
    boot.to_csv(args.output_dir / "jshis_official_hazard_response_bootstrap_ci.csv", index=False)
    write_report(
        frame,
        fold_summary,
        boot,
        args.output_dir / "jshis_official_hazard_response_uncertainty_spatial_check.md",
    )
    build_figure(frame, fold_summary, boot, args.figure_dir, args.figure_stem, args.language)

    print(f"Wrote {args.output_dir / 'jshis_official_hazard_response_fold_summary.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_official_hazard_response_bootstrap_ci.csv'}")
    print(f"Wrote {args.figure_dir / (args.figure_stem + '.pdf')}")


if __name__ == "__main__":
    main()
