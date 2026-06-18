#!/usr/bin/env python3
"""Build the main-text spatial evidence figure for the CEE manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    HAS_CARTOPY = True
except Exception:
    HAS_CARTOPY = False


HERE = Path(__file__).resolve().parent
SUPPLEMENT = HERE / "supplement"
OUTDIR = HERE / "figures"

SURFACE_MESH = SUPPLEMENT / "jshis_sa3_continuous_correction_surface_mesh.csv"
SURFACE_SUMMARY = SUPPLEMENT / "jshis_sa3_continuous_correction_surface_summary.csv"
SURFACE_STATIONS = SUPPLEMENT / "jshis_sa3_continuous_correction_surface_station_validation.csv"
VS400_MESH = SUPPLEMENT / "jshis_zamp_vs400_site_amplification_mesh.csv"
VS400_SUMMARY = SUPPLEMENT / "jshis_zamp_vs400_site_amplification_summary.csv"
OFFICIAL_SA3 = SUPPLEMENT / "jshis_official_response_sa3_station_values.csv"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "font.size": 8.8,
            "axes.titlesize": 9.6,
            "axes.labelsize": 8.8,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.8,
            "legend.fontsize": 7.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 300,
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
        fontsize=11.4,
        fontweight="bold",
        color="#111111",
    )
    ax.text(
        -0.005,
        1.04,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.2,
        color="#111111",
    )


def map_axis(fig: plt.Figure, spec):
    if HAS_CARTOPY:
        ax = fig.add_subplot(spec, projection=ccrs.PlateCarree())
        ax.set_extent([122.0, 149.5, 23.5, 46.8], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#F6F6F3", edgecolor="none", zorder=0)
        ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#F8FBFC", edgecolor="none", zorder=0)
        ax.coastlines(resolution="50m", linewidth=0.45, color="#666666", zorder=2)
        gl = ax.gridlines(draw_labels=True, linewidth=0.25, color="#D2D2D2", alpha=0.8)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": 7.2}
        gl.ylabel_style = {"size": 7.2}
        return ax

    ax = fig.add_subplot(spec)
    ax.set_xlim(122.0, 149.5)
    ax.set_ylim(23.5, 46.8)
    ax.grid(color="#E6E6E6", lw=0.35)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    return ax


def scatter_map(ax, df: pd.DataFrame, values: pd.Series, cmap, norm, size: float = 7.0):
    kwargs = dict(
        x=df["lon"],
        y=df["lat"],
        c=values,
        cmap=cmap,
        norm=norm,
        s=size,
        lw=0,
        alpha=0.92,
        rasterized=True,
        zorder=3,
    )
    if HAS_CARTOPY:
        kwargs["transform"] = ccrs.PlateCarree()
    return ax.scatter(**kwargs)


def add_colorbar(fig: plt.Figure, ax, mappable, label: str) -> None:
    cbar = fig.colorbar(mappable, ax=ax, fraction=0.038, pad=0.02, shrink=0.86)
    cbar.set_label(label, fontsize=7.6)
    cbar.ax.tick_params(labelsize=7.1, length=2)


def group_row(summary: pd.DataFrame, group: str) -> pd.Series:
    row = summary.loc[summary["group"].eq(group)]
    if row.empty:
        raise ValueError(f"Missing group {group}")
    return row.iloc[0]


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    setup_style()

    surface = pd.read_csv(SURFACE_MESH)
    surface_summary = pd.read_csv(SURFACE_SUMMARY)
    stations = pd.read_csv(SURFACE_STATIONS)
    vs400 = pd.read_csv(VS400_MESH)
    vs400_summary = pd.read_csv(VS400_SUMMARY)
    official = pd.read_csv(OFFICIAL_SA3)

    q10 = official.loc[official["probability_level"].eq("50y_10pct")].copy()

    overall = surface_summary.loc[
        surface_summary["row_type"].eq("overall_station_validation_site_space_hgb")
    ].iloc[0]
    surface_groups = surface_summary.loc[surface_summary["row_type"].str.startswith("mesh_group_", na=False)].copy()
    surface_groups["group_label"] = surface_groups["row_type"].str.replace("mesh_group_", "", regex=False)

    fig = plt.figure(figsize=(7.35, 7.85), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.12, 1.02, 0.98])

    ax_a = map_axis(fig, gs[0, 0])
    norm_a = TwoSlopeNorm(vmin=-0.32, vcenter=0.0, vmax=0.08)
    sc_a = scatter_map(
        ax_a,
        surface,
        surface["surface_sa3_correction_log10"],
        cmap="RdBu_r",
        norm=norm_a,
        size=7.2,
    )
    tail = surface.loc[surface["is_tail"].astype(bool)]
    if not tail.empty:
        map_kwargs = dict(x=tail["lon"], y=tail["lat"], s=5.2, facecolors="none", edgecolors="#222222", lw=0.20, alpha=0.36, zorder=4)
        if HAS_CARTOPY:
            map_kwargs["transform"] = ccrs.PlateCarree()
        ax_a.scatter(**map_kwargs)
    panel_label(ax_a, "A", "SA(3.0 s) sampled-grid correction")
    add_colorbar(fig, ax_a, sc_a, "log10 correction")

    ax_b = map_axis(fig, gs[0, 1])
    norm_b = Normalize(vmin=np.nanpercentile(vs400["zamp_avs30_mps"], 2), vmax=np.nanpercentile(vs400["zamp_avs30_mps"], 98))
    sc_b = scatter_map(ax_b, vs400, vs400["zamp_avs30_mps"], cmap="viridis", norm=norm_b, size=7.0)
    tail_vs = vs400.loc[vs400["is_tail"].astype(bool)]
    if not tail_vs.empty:
        low_kwargs = dict(
            x=tail_vs["lon"],
            y=tail_vs["lat"],
            s=5.0,
            facecolors="none",
            edgecolors="#222222",
            alpha=0.30,
            lw=0.20,
            zorder=4,
        )
        if HAS_CARTOPY:
            low_kwargs["transform"] = ccrs.PlateCarree()
        ax_b.scatter(**low_kwargs)
    panel_label(ax_b, "B", "J-SHIS VS400 gridded AVS30")
    add_colorbar(fig, ax_b, sc_b, "AVS30 (m/s)")

    ax_c = fig.add_subplot(gs[1, 0])
    station_size = 8 + 42 * np.sqrt(np.clip(stations["n_records"], 1, None)) / np.sqrt(stations["n_records"].max())
    ax_c.scatter(
        stations["observed_sa3_correction_log10"],
        stations["site_space_predicted_sa3_correction_log10"],
        s=station_size,
        color="#4C78A8",
        alpha=0.32,
        edgecolors="none",
        rasterized=True,
    )
    lo = float(np.nanmin([stations["observed_sa3_correction_log10"].min(), stations["site_space_predicted_sa3_correction_log10"].min()]))
    hi = float(np.nanmax([stations["observed_sa3_correction_log10"].max(), stations["site_space_predicted_sa3_correction_log10"].max()]))
    pad = 0.05
    ax_c.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#111111", lw=0.8)
    ax_c.set_xlim(lo - pad, hi + pad)
    ax_c.set_ylim(lo - pad, hi + pad)
    ax_c.set_xlabel("Held-out station correction (log10)")
    ax_c.set_ylabel("Predicted correction (log10)")
    ax_c.grid(color="#E6E6E6", lw=0.45)
    panel_label(ax_c, "C", "Spatial-block station validation")
    ax_c.text(
        0.03,
        0.97,
        f"RMSE {overall['weighted_rmse']:.3f}\n"
        f"vs zero {overall['rmse_reduction_vs_zero_pct']:.1f}%\n"
        f"vs national mean {overall['rmse_reduction_vs_mean_pct']:.1f}%",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#DDDDDD", "alpha": 0.88},
    )

    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.scatter(
        q10["official_sa3_g"],
        q10["corrected_sa3_g"],
        s=10,
        color="#4C78A8",
        alpha=0.42,
        edgecolors="none",
        rasterized=True,
    )
    limit = float(np.nanpercentile(q10[["official_sa3_g", "corrected_sa3_g"]].to_numpy(), 99.8)) * 1.1
    ax_d.plot([0, limit], [0, limit], color="#111111", lw=0.8)
    ax_d.set_xlim(0, limit)
    ax_d.set_ylim(0, limit)
    ax_d.set_xlabel("Official SA(3.0 s), 50-year 10% (g)")
    ax_d.set_ylabel("Station-corrected SA(3.0 s) (g)")
    ax_d.grid(color="#E6E6E6", lw=0.45)
    panel_label(ax_d, "D", "Official mesh ordinate sensitivity")
    ax_d.text(
        0.03,
        0.97,
        "median 0.076 to 0.048 g",
        transform=ax_d.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#DDDDDD", "alpha": 0.88},
    )

    ax_e = fig.add_subplot(gs[2, 0])
    bins = pd.qcut(vs400["zamp_avs30_mps"], 4, labels=["Q1 soft", "Q2", "Q3", "Q4 stiff"], duplicates="drop")
    tail_by_bin = vs400.assign(avs_bin=bins).groupby("avs_bin", observed=True)["is_tail"].mean().mul(100)
    bars = ax_e.bar(tail_by_bin.index.astype(str), tail_by_bin.values, color=["#B23A48", "#6B8FB6", "#6B8FB6", "#6B8FB6"], width=0.68)
    ax_e.bar_label(bars, labels=[f"{v:.1f}%" for v in tail_by_bin.values], padding=2, fontsize=7.6)
    ax_e.set_ylabel("Residual-tail meshes (%)")
    ax_e.set_ylim(0, max(tail_by_bin.values) * 1.28)
    ax_e.grid(axis="y", color="#E6E6E6", lw=0.45)
    panel_label(ax_e, "E", "Residual tails by gridded AVS30")

    ax_f = fig.add_subplot(gs[2, 1])
    selected = [
        ("All", "all", "all_meshes"),
        ("Tail", "tail", "residual_tail"),
        ("NE Hokkaido", "ne_hokkaido", "northeast_hokkaido"),
        ("D1400 1-3 s", "d1400_1_3s_window", "d1400_1_3s_window"),
    ]
    y = np.arange(len(selected))
    medians = np.array([group_row(vs400_summary, vs_key)["surface_multiplier_median"] for _, vs_key, _ in selected])
    p05 = np.array(
        [
            surface_groups.loc[surface_groups["group_label"].eq(surface_key), "surface_multiplier_p05"].iloc[0]
            for _, _, surface_key in selected
        ]
    )
    p95 = np.array(
        [
            surface_groups.loc[surface_groups["group_label"].eq(surface_key), "surface_multiplier_p95"].iloc[0]
            for _, _, surface_key in selected
        ]
    )
    ax_f.hlines(y, p05, p95, color="#9A9A9A", lw=1.9)
    ax_f.plot(medians, y, "o", color="#4C78A8", ms=4.5)
    ax_f.axvline(1.0, color="#222222", lw=0.8, ls="--")
    ax_f.set_yticks(y, [label for label, _, _ in selected])
    ax_f.set_xlim(0.42, 1.03)
    ax_f.set_xlabel("SA(3.0 s) multiplier")
    ax_f.grid(axis="x", color="#E6E6E6", lw=0.45)
    panel_label(ax_f, "F", "Sampled-grid multiplier ranges")
    ax_f.legend(
        handles=[
            Line2D([0], [0], color="#9A9A9A", lw=1.9, label="5%-95%"),
            Line2D([0], [0], marker="o", color="#4C78A8", lw=0, label="median"),
        ],
        frameon=False,
        loc="lower right",
    )

    for ax in [ax_c, ax_d, ax_e, ax_f]:
        ax.tick_params(length=3)

    fig.savefig(OUTDIR / "figure5_spatial_evidence.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure5_spatial_evidence.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
