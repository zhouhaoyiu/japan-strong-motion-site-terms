#!/usr/bin/env python3
"""Build a publication-style overview figure for the English manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "figures"
PROJECT_ROOT = HERE.parents[1]

SITE_SCHEMA = PROJECT_ROOT / "outputs" / "jshis_site_schema_v2024_sub1.csv"
SOURCE_SCHEMA = PROJECT_ROOT / "outputs" / "jshis_source_schema_v2024_sub1.csv"
NETWORK_COUNTS = PROJECT_ROOT / "outputs" / "jshis_smrec_schema_sub1_network_counts.csv"
SITE_EFFECTS = PROJECT_ROOT / "outputs" / "jshis_mf2013_official_site_correction_effects.csv"
STATION_MODEL = PROJECT_ROOT / "outputs" / "jshis_station_model_ablation_summary.csv"
HAZARD_SUMMARY = PROJECT_ROOT / "outputs" / "jshis_official_hazard_response_summary.csv"
SIGMA_METRICS = PROJECT_ROOT / "outputs" / "jshis_public_psha_sigma_sensitivity1797_metrics.csv"
EXECUTION_BY_PROBABILITY = PROJECT_ROOT / "outputs" / "jshis_public_psha_national_stratified10467_arrays8_recompute_by_probability.csv"
EXECUTION_BY_GEO = PROJECT_ROOT / "outputs" / "jshis_public_psha_national_stratified10467_arrays8_recompute_by_coarse_geo.csv"

BLUE = "#4C78A8"
ORANGE = "#D08A39"
TEAL = "#4F9A8A"
RED = "#C44E52"
GRAY = "#4A4A4A"
TEXT = "#202020"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "font.size": 8.8,
            "axes.labelsize": 8.8,
            "axes.titlesize": 9.0,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.8,
            "legend.fontsize": 7.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def panel_label(ax, label: str, title: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.5,
        fontweight="bold",
        color="#111111",
    )
    ax.text(
        0.01,
        1.04,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.9,
        fontweight="bold",
        color="#111111",
    )


def draw_map(ax, sites: pd.DataFrame, sources: pd.DataFrame) -> None:
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        ax.set_extent([122.5, 149.0, 23.5, 46.5], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#F2F2F0", edgecolor="none", zorder=0)
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.45, edgecolor="#777777", zorder=1)
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.25, edgecolor="#A0A0A0", zorder=1)
        transform = ccrs.PlateCarree()
        plot_kwargs = {"transform": transform}
        gridlines = ax.gridlines(
            draw_labels=True,
            linewidth=0.28,
            color="#D8D8D8",
            linestyle="-",
            x_inline=False,
            y_inline=False,
        )
        gridlines.top_labels = False
        gridlines.right_labels = False
        gridlines.xlabel_style = {"size": 7}
        gridlines.ylabel_style = {"size": 7}
    except Exception:
        ax.set_xlim(122.5, 149.0)
        ax.set_ylim(23.5, 46.5)
        ax.set_xlabel("Longitude (deg E)")
        ax.set_ylabel("Latitude (deg N)")
        ax.grid(color="#E6E6E6", lw=0.45)
        plot_kwargs = {}

    source_plot = sources.dropna(subset=["jem_lon", "jem_lat", "mjma"]).copy()
    source_plot = source_plot[source_plot["mjma"] >= 5.0]
    sizes = np.clip((source_plot["mjma"].astype(float) - 4.7) ** 2 * 4.5, 3.0, 80.0)
    ax.scatter(
        source_plot["jem_lon"],
        source_plot["jem_lat"],
        s=sizes,
        color="#B7B7B7",
        alpha=0.38,
        edgecolors="none",
        rasterized=True,
        zorder=2,
        **plot_kwargs,
    )
    for network_id, color in [(1, BLUE), (2, ORANGE)]:
        sub = sites[sites["obs_network_id"].eq(network_id)]
        ax.scatter(
            sub["lon"],
            sub["lat"],
            s=5.8,
            color=color,
            alpha=0.72,
            edgecolors="white",
            linewidths=0.10,
            rasterized=True,
            zorder=3,
            **plot_kwargs,
        )
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white", markersize=4.5, label="K-NET station"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=ORANGE, markeredgecolor="white", markersize=4.5, label="KiK-net station"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#B7B7B7", markeredgecolor="none", markersize=4.5, label=r"$M\geq5$ source"),
        ],
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.84,
        loc="lower left",
        borderpad=0.25,
        labelspacing=0.25,
        handletextpad=0.35,
    )
    panel_label(ax, "A", "Record and source coverage")


def draw_design(ax, network_counts: pd.DataFrame) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    panel_label(ax, "B", "Evidence chain")

    total_records = int(network_counts["n_smrec_records"].sum())
    knet = int(network_counts.loc[network_counts["network_label"].eq("K-NET"), "n_smrec_records"].iloc[0])
    kik = int(network_counts.loc[network_counts["network_label"].eq("KiK-net"), "n_smrec_records"].iloc[0])
    total = knet + kik

    ax.text(0.03, 0.88, f"{total_records:,}", ha="left", va="center", fontsize=18.5, fontweight="bold", color=TEXT)
    ax.text(0.03, 0.77, "public strong-motion records", ha="left", va="center", fontsize=8.4, color=TEXT)
    ax.text(0.03, 0.67, "2,581 stations   1,840 sources", ha="left", va="center", fontsize=7.9, color=GRAY)

    y0 = 0.56
    ax.add_patch(Rectangle((0.03, y0), 0.54, 0.032, facecolor="#F0F0F0", edgecolor="none"))
    ax.add_patch(Rectangle((0.03, y0), 0.54 * knet / total, 0.032, facecolor=BLUE, edgecolor="none"))
    ax.add_patch(Rectangle((0.03 + 0.54 * knet / total, y0), 0.54 * kik / total, 0.032, facecolor=ORANGE, edgecolor="none"))
    ax.text(0.03, y0 - 0.041, "K-NET 38.6%", fontsize=7.1, ha="left", va="center", color=BLUE)
    ax.text(0.57, y0 - 0.041, "KiK-net 61.4%", fontsize=7.1, ha="right", va="center", color=ORANGE)

    rows = [
        ("Residuals", "MF2013 basic + site terms"),
        ("Holdout", "spatial-block validation"),
        ("Spectrum", "J-SHIS SA(3.0 s) ordinates"),
        ("PSHA audit", "10,467 meshes; sigma; tails"),
    ]
    x_left, x_mid, x_right = 0.03, 0.45, 0.96
    y_top, row_h = 0.43, 0.115
    ax.plot([x_left, x_right], [y_top + 0.035, y_top + 0.035], color="#1F1F1F", lw=0.8)
    ax.plot([x_left, x_right], [y_top - len(rows) * row_h - 0.012, y_top - len(rows) * row_h - 0.012], color="#1F1F1F", lw=0.8)
    for i, (label, value) in enumerate(rows):
        y = y_top - i * row_h
        ax.text(x_left, y, label, ha="left", va="center", fontsize=7.2, fontweight="bold", color=TEXT)
        ax.text(x_mid, y, value, ha="left", va="center", fontsize=7.0, color=TEXT)
        if i < len(rows) - 1:
            ax.plot([x_left, x_right], [y - row_h / 2, y - row_h / 2], color="#E1E1E1", lw=0.6)


def draw_period_gain(ax, site_effects: pd.DataFrame, station_model: pd.DataFrame) -> None:
    order = ["PGA RotD50", "SA(0.3s) RotD50", "SA(1.0s) RotD50", "SA(3.0s) RotD50"]
    labels = ["PGA", "0.3", "1.0", "3.0"]
    effect = site_effects[
        site_effects["effect_type"].eq("mae")
        & site_effects["feature"].eq("overall")
        & site_effects["target_label"].isin(order)
    ].copy()
    effect["target_label"] = pd.Categorical(effect["target_label"], categories=order, ordered=True)
    effect = effect.sort_values("target_label")

    model = station_model[
        station_model["split"].eq("spatial_block")
        & station_model["feature_set"].eq("full_site_location")
        & station_model["target_label"].isin(order)
    ].copy()
    model["target_label"] = pd.Categorical(model["target_label"], categories=order, ordered=True)
    model = model.sort_values("target_label")

    x = np.arange(len(order))
    width = 0.34
    ax.bar(x - width / 2, effect["relative_reduction_pct"].astype(float), width=width, color=BLUE, label="MF2013 site term")
    ax.bar(x + width / 2, model["rmse_reduction_pct_mean"].astype(float), width=width, color=TEAL, label="held-out station model")
    ax.axhline(0, color="#333333", lw=0.75)
    ax.set_xticks(x, labels)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Error reduction (%)")
    ax.set_ylim(-6, 52)
    ax.grid(axis="y", color="#E6E6E6", lw=0.6)
    ax.legend(frameon=False, loc="upper left", handlelength=1.5, labelspacing=0.25)
    panel_label(ax, "C", "Period-dependent reduction")


def draw_hazard(ax, hazard_summary: pd.DataFrame) -> None:
    row = hazard_summary[hazard_summary["metric"].eq("official_sa3_50y_10pct")].iloc[0]
    quantiles = ["5%", "50%", "95%"]
    official = np.array([row["official_sa3_g_q05"], row["official_sa3_g_q50"], row["official_sa3_g_q95"]], dtype=float)
    corrected = np.array([row["corrected_sa3_g_q05"], row["corrected_sa3_g_q50"], row["corrected_sa3_g_q95"]], dtype=float)
    x = np.arange(len(quantiles))
    ax.plot(x, official, color="#333333", lw=1.4, marker="o", ms=4.2, label="official ordinate")
    ax.plot(x, corrected, color=RED, lw=1.4, marker="o", ms=4.2, label="station-corrected")
    for xi, y0, y1 in zip(x, official, corrected, strict=True):
        ax.plot([xi, xi], [y1, y0], color="#BDBDBD", lw=0.9, zorder=0)
    ax.set_xticks(x, quantiles)
    ax.set_xlabel("Station-distribution percentile")
    ax.set_ylabel("SA(3.0 s), 50-year 10% (g)")
    ax.set_ylim(0, max(official.max(), corrected.max()) * 1.22)
    ax.grid(axis="y", color="#E6E6E6", lw=0.6)
    ax.legend(frameon=False, loc="upper left", handlelength=1.5, labelspacing=0.25)
    ax.text(
        1,
        corrected[1] * 0.70,
        f"median {official[1]:.3f} to {corrected[1]:.3f} g",
        ha="center",
        va="top",
        fontsize=7.5,
        color=TEXT,
    )
    panel_label(ax, "D", "Official SA(3.0 s) ordinates")


def draw_sigma(ax, sigma_metrics: pd.DataFrame) -> None:
    all_rows = sigma_metrics[sigma_metrics["probability_column"].eq("ALL")].sort_values("sigma_log10")
    ax.plot(
        all_rows["sigma_log10"],
        all_rows["rmse"],
        marker="o",
        color=GRAY,
        lw=1.7,
        label="all thresholds",
    )
    for col, color in [("T50_P02_BV", BLUE), ("T50_P05_BV", ORANGE), ("T50_P10_BV", TEAL), ("T50_P39_BV", RED)]:
        sub = sigma_metrics[sigma_metrics["probability_column"].eq(col)].sort_values("sigma_log10")
        ax.plot(sub["sigma_log10"], sub["rmse"], marker="o", ms=3.0, lw=0.95, color=color, alpha=0.68)
    best = all_rows.loc[all_rows["rmse"].idxmin()]
    ax.axvline(float(best["sigma_log10"]), color="#111111", lw=0.8, ls="--")
    ax.text(
        float(best["sigma_log10"]) + 0.006,
        float(best["rmse"]) + 0.010,
        f"global RMSE min\nsigma={best['sigma_log10']:.2f}",
        ha="left",
        va="bottom",
        fontsize=7.1,
        color=TEXT,
    )
    ax.set_xlabel(r"log$_{10}$ sigma")
    ax.set_ylabel("RMSE of log10 residual")
    ax.set_ylim(0.18, max(sigma_metrics["rmse"]) * 1.13)
    ax.grid(axis="y", color="#E6E6E6", lw=0.6)
    panel_label(ax, "E", "Sigma audit")


def draw_root_cause(ax, by_probability: pd.DataFrame, by_geo: pd.DataFrame) -> None:
    p = by_probability.sort_values("probability_column").copy()
    labels = ["2%", "5%", "10%", "39%"]
    ax.bar(np.arange(len(p)), p["rmse"], color=[BLUE, ORANGE, TEAL, RED], width=0.62)
    ax.set_xticks(np.arange(len(p)), labels)
    ax.set_xlabel("50-year exceedance probability")
    ax.set_ylabel("RMSE of log10 residual")
    ax.grid(axis="y", color="#E6E6E6", lw=0.6)
    ax.set_ylim(0, max(p["rmse"]) * 1.34)
    worst_geo = by_geo.sort_values("abs_p95", ascending=False).iloc[0]
    ax.text(
        0.02,
        0.97,
        "10,467 meshes\n41,868 valid rows\nzero-probability rows: 0",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        color=TEXT,
        bbox={"boxstyle": "round,pad=0.24", "facecolor": "white", "edgecolor": "#D9D9D9", "linewidth": 0.6},
    )
    geo_label = {
        "northeast_hokkaido_band": "NE Hokkaido",
        "low_latitude_islands": "low-latitude islands",
        "southwest_japan_band": "SW Japan",
        "central_japan_band": "central Japan",
        "tohoku_hokkaido_band": "Tohoku-Hokkaido",
    }.get(str(worst_geo["coarse_geo_band"]), str(worst_geo["coarse_geo_band"]).replace("_", " "))
    ax.text(
        0.98,
        0.74,
        f"largest diagnostic tail\n{geo_label}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
        color=TEXT,
    )
    panel_label(ax, "F", "Residual-tail audit")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    setup_style()
    sites = pd.read_csv(SITE_SCHEMA)
    sources = pd.read_csv(SOURCE_SCHEMA)
    network_counts = pd.read_csv(NETWORK_COUNTS)
    site_effects = pd.read_csv(SITE_EFFECTS)
    station_model = pd.read_csv(STATION_MODEL)
    hazard_summary = pd.read_csv(HAZARD_SUMMARY)
    sigma_metrics = pd.read_csv(SIGMA_METRICS)
    root_by_probability = pd.read_csv(EXECUTION_BY_PROBABILITY)
    root_by_geo = pd.read_csv(EXECUTION_BY_GEO)

    try:
        import cartopy.crs as ccrs

        fig = plt.figure(figsize=(12.8, 7.0), dpi=300)
        grid = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.0, 1.0], height_ratios=[1.03, 0.97], hspace=0.43, wspace=0.34)
        ax_a = fig.add_subplot(grid[0, 0], projection=ccrs.PlateCarree())
    except Exception:
        fig = plt.figure(figsize=(12.8, 7.0), dpi=300)
        grid = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.0, 1.0], height_ratios=[1.03, 0.97], hspace=0.43, wspace=0.34)
        ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[1, 0])
    ax_e = fig.add_subplot(grid[1, 1])
    ax_f = fig.add_subplot(grid[1, 2])

    draw_map(ax_a, sites, sources)
    draw_design(ax_b, network_counts)
    draw_period_gain(ax_c, site_effects, station_model)
    draw_hazard(ax_d, hazard_summary)
    draw_sigma(ax_e, sigma_metrics)
    draw_root_cause(ax_f, root_by_probability, root_by_geo)

    fig.savefig(OUTDIR / "figure1.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure1.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
