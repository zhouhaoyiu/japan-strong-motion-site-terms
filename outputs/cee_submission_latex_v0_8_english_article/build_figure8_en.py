#!/usr/bin/env python3
"""Build English Figure 8 from official J-SHIS and KiK-net audit outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, FixedLocator
from matplotlib.transforms import Bbox


HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "figures"
PROJECT_ROOT = HERE.parents[1]
STATION_IMPACT = PROJECT_ROOT / "outputs" / "jshis_station_hazard_impact_by_station.csv"
RESPONSE_SUMMARY = PROJECT_ROOT / "outputs" / "jshis_official_hazard_response_summary.csv"
SA3_VALUES = PROJECT_ROOT / "outputs" / "jshis_official_response_sa3_station_values.csv"
UHS_EXAMPLES = PROJECT_ROOT / "outputs" / "jshis_official_response_uhs_examples.csv"
KIKNET_RATIOS = PROJECT_ROOT / "outputs" / "kiknet_multievent_transfer_functions.csv"
PHYSICAL_LINK = PROJECT_ROOT / "outputs" / "jshis_public_psha_10467_near_station_physical_link.csv"
PHYSICAL_SUMMARY = PROJECT_ROOT / "outputs" / "jshis_public_psha_10467_near_station_physical_summary.csv"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "font.size": 10.8,
            "axes.titlesize": 12.0,
            "axes.labelsize": 10.6,
            "xtick.labelsize": 9.2,
            "ytick.labelsize": 9.2,
            "legend.fontsize": 8.4,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def label_for(example_type: str) -> str:
    return {
        "maximum_negative": "maximum negative",
        "q50": "median",
        "maximum_positive": "maximum positive",
    }.get(example_type, example_type)


def log_plain(value: float, _position: int) -> str:
    if value <= 0:
        return ""
    exponent = int(np.round(np.log10(value)))
    if np.isclose(value, 10.0**exponent):
        return f"1e{exponent}"
    return ""


def probability_label(metric: str) -> str:
    label = metric.replace("official_sa3_50y_", "").replace("pct", "")
    return label


def select_summary(summary: pd.DataFrame) -> pd.DataFrame:
    rows = summary[summary["metric"].astype(str).str.startswith("official_sa3_50y_")].copy()
    rows["probability_pct"] = rows["metric"].map(lambda value: float(probability_label(str(value))))
    return rows.sort_values("probability_pct")


def panel_title(ax, label: str, title: str) -> None:
    ax.text(
        -0.10,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=13.2,
        fontweight="bold",
        color="#111111",
    )
    ax.text(
        -0.035,
        1.04,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.4,
        color="#111111",
    )


def save_axes_subset(fig, axes, path: Path) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bbox = Bbox.union([ax.get_tightbbox(renderer) for ax in axes])
    bbox = bbox.expanded(1.04, 1.06).transformed(fig.dpi_scale_trans.inverted())
    kwargs = {"dpi": 240} if path.suffix.lower() == ".png" else {}
    fig.savefig(path, bbox_inches=bbox, **kwargs)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    setup_style()
    station_impact = pd.read_csv(STATION_IMPACT)
    response_summary = pd.read_csv(RESPONSE_SUMMARY)
    sa3_values = pd.read_csv(SA3_VALUES)
    uhs_examples = pd.read_csv(UHS_EXAMPLES)
    kiknet = pd.read_csv(KIKNET_RATIOS)
    physical_link = pd.read_csv(PHYSICAL_LINK)
    physical_summary = pd.read_csv(PHYSICAL_SUMMARY)
    uhs_examples["period_s"] = pd.to_numeric(uhs_examples["period_s"], errors="coerce")
    uhs_examples["official_sa_g"] = pd.to_numeric(uhs_examples["official_sa_g"], errors="coerce")
    uhs_examples["station_corrected_sa_g"] = pd.to_numeric(uhs_examples["station_corrected_sa_g"], errors="coerce")

    fig = plt.figure(figsize=(12.8, 11.2))
    grid = fig.add_gridspec(3, 2, hspace=0.58, wspace=0.34)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    ax_e = fig.add_subplot(grid[2, 0])
    ax_f = fig.add_subplot(grid[2, 1])

    factors = station_impact["sa3_uhs_multiplier"].astype(float)
    ax_a.hist(
        factors,
        bins=np.geomspace(max(factors.min(), 0.08), factors.max() * 1.05, 34),
        color="#4C78A8",
        edgecolor="white",
    )
    ax_a.set_xscale("log")
    ax_a.set_xlim(0.1, 1.8)
    ax_a.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1.0, 1.5]))
    ax_a.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax_a.axvline(1.0, color="black", lw=1.0)
    for q, ls in [(0.05, "--"), (0.50, "-"), (0.95, "--")]:
        ax_a.axvline(factors.quantile(q), color="#C44E52", lw=1.0, ls=ls)
    ax_a.set_xlabel("Station multiplier at SA(3.0 s)")
    ax_a.set_ylabel("Number of stations")
    panel_title(ax_a, "A", "Station multipliers from held-out residuals")
    ax_a.grid(axis="y", color="#E6E6E6")
    ax_a.legend(
        handles=[
            Line2D([0], [0], color="black", lw=1.0, label="multiplier = 1"),
            Line2D([0], [0], color="#C44E52", lw=1.0, ls="-", label="median"),
            Line2D([0], [0], color="#C44E52", lw=1.0, ls="--", label="5% / 95%"),
        ],
        frameon=False,
        loc="upper left",
    )

    colors = {"maximum_negative": "#4C78A8", "q50": "#555555", "maximum_positive": "#C44E52"}
    summary = select_summary(response_summary)
    x = summary["probability_pct"].to_numpy()
    official_q05 = summary["official_sa3_g_q05"].to_numpy()
    official_q50 = summary["official_sa3_g_q50"].to_numpy()
    official_q95 = summary["official_sa3_g_q95"].to_numpy()
    corrected_q05 = summary["corrected_sa3_g_q05"].to_numpy()
    corrected_q50 = summary["corrected_sa3_g_q50"].to_numpy()
    corrected_q95 = summary["corrected_sa3_g_q95"].to_numpy()
    ax_b.fill_between(x, official_q05, official_q95, color="#8C8C8C", alpha=0.16, lw=0)
    ax_b.fill_between(x, corrected_q05, corrected_q95, color="#4C78A8", alpha=0.16, lw=0)
    ax_b.plot(x, official_q50, color="#333333", lw=1.9, marker="o", label="official median")
    ax_b.plot(x, corrected_q50, color="#4C78A8", lw=1.9, marker="o", label="corrected median")
    ax_b.set_xticks(x)
    ax_b.set_xlabel("50-year exceedance probability (%)")
    ax_b.set_ylabel("SA(3.0 s) ordinate (g)")
    panel_title(ax_b, "B", "Official spectrum ordinate quantiles")
    ax_b.grid(color="#E6E6E6")
    ax_b.legend(
        handles=[
            Line2D([0], [0], color="#333333", lw=1.9, marker="o", label="official median"),
            Line2D([0], [0], color="#4C78A8", lw=1.9, marker="o", label="corrected median"),
            Line2D([0], [0], color="#8C8C8C", lw=7, alpha=0.20, label="5%-95% range"),
        ],
        frameon=False,
        loc="upper right",
    )

    q10 = sa3_values[sa3_values["probability_level"].eq("50y_10pct")].copy()
    ax_c.scatter(
        q10["official_sa3_g"],
        q10["corrected_sa3_g"],
        s=10,
        color="#4C78A8",
        alpha=0.42,
        edgecolors="none",
    )
    limit = float(np.nanmax([q10["official_sa3_g"].max(), q10["corrected_sa3_g"].max()])) * 1.05
    ax_c.plot([0, limit], [0, limit], color="black", lw=0.9)
    ax_c.set_xlim(0, limit)
    ax_c.set_ylim(0, limit)
    ax_c.set_xlabel("Official SA(3.0 s), 50-year 10% (g)")
    ax_c.set_ylabel("Station-corrected SA(3.0 s) (g)")
    panel_title(ax_c, "C", "Matched station-mesh ordinates")
    ax_c.grid(color="#E6E6E6")

    uhs = uhs_examples[
        uhs_examples["probability_level"].eq("50y_10pct")
        & uhs_examples["example_type"].astype(str).isin(["maximum_negative", "q50", "maximum_positive"])
    ].copy()
    official_styles = {
        "maximum_negative": (0, (5, 2)),
        "q50": (0, (3, 2)),
        "maximum_positive": (0, (1.4, 1.8)),
    }
    for example_type in ["maximum_negative", "q50", "maximum_positive"]:
        sub = uhs[uhs["example_type"].astype(str).eq(example_type)].sort_values("period_s")
        if sub.empty:
            continue
        marker_index = int(np.argmin(np.abs(sub["period_s"].to_numpy() - 3.0)))
        ax_d.plot(
            sub["period_s"],
            sub["station_corrected_sa_g"],
            color=colors[example_type],
            lw=1.7,
            marker="o",
            ms=3.0,
            markevery=[marker_index],
            label=f"{label_for(example_type)}: {sub['site_code'].iloc[0]}",
            zorder=2,
        )
        ax_d.plot(
            sub["period_s"],
            sub["official_sa_g"],
            color="#222222",
            lw=1.05,
            ls=official_styles[example_type],
            alpha=0.82,
            marker="o",
            ms=3.0,
            mfc="white",
            mec="#222222",
            markevery=[marker_index],
            zorder=3,
        )
    ax_d.set_xscale("log")
    ax_d.set_xlim(0.08, 6.0)
    ax_d.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0]))
    ax_d.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax_d.set_xlabel("Period (s)")
    ax_d.set_ylabel("5%-damped SA, 50-year 10% (g)")
    panel_title(ax_d, "D", "Uniform-hazard spectra after SA(3.0 s) shift")
    ax_d.axvline(3.0, color="#777777", lw=0.75, ls=":", zorder=1)
    ax_d.grid(color="#E6E6E6", which="both")
    station_legend = ax_d.legend(frameon=False, loc="upper right")
    ax_d.add_artist(station_legend)
    ax_d.legend(
        handles=[
            Line2D([0], [0], color="#222222", lw=1.05, ls="--", label="official spectrum"),
            Line2D([0], [0], color="#333333", lw=1.7, ls="-", label="SA(3.0 s) adjusted"),
        ],
        frameon=False,
        loc="lower left",
    )

    band_specs = [
        ("0.2-0.5", "spectral_ratio_0.2_0.5_hz"),
        ("0.5-1", "spectral_ratio_0.5_1_hz"),
        ("1-2", "spectral_ratio_1_2_hz"),
        ("2-5", "spectral_ratio_2_5_hz"),
        ("5-10", "spectral_ratio_5_10_hz"),
        ("10-20", "spectral_ratio_10_20_hz"),
    ]
    box_data = []
    box_labels = []
    medians = []
    for label, col in band_specs:
        vals = pd.to_numeric(kiknet[col], errors="coerce")
        vals = vals[np.isfinite(vals) & (vals > 0)]
        box_data.append(vals.to_numpy())
        box_labels.append(label)
        medians.append(float(vals.median()))
    box = ax_e.boxplot(box_data, showfliers=False, patch_artist=True, widths=0.58)
    for patch in box["boxes"]:
        patch.set(facecolor="#88B0D8", edgecolor="#355F8A", linewidth=0.9)
    for item in box["medians"]:
        item.set(color="#111111", linewidth=1.2)
    for item in box["whiskers"] + box["caps"]:
        item.set(color="#355F8A", linewidth=0.9)
    ax_e.axhline(1.0, color="#111111", lw=0.95)
    ax_e.set_yscale("log")
    ax_e.set_ylim(0.18, 5.8)
    ax_e.set_xticks(np.arange(1, len(box_labels) + 1), box_labels)
    ax_e.set_xlabel("Frequency band (Hz)")
    ax_e.set_ylabel("Surface/downhole spectral ratio")
    panel_title(ax_e, "E", "Multi-event KiK-net borehole/surface check")
    ax_e.grid(axis="y", color="#E6E6E6", which="both")
    for i, med in enumerate(medians, start=1):
        ax_e.text(i, med * 1.06, f"{med:.2g}", ha="center", va="bottom", fontsize=7.6)

    all_meshes = physical_link.copy()
    tail = all_meshes[all_meshes["is_tail"].astype(bool)]
    ne_tail = all_meshes[all_meshes["is_ne_hokkaido"].astype(bool) & all_meshes["is_tail"].astype(bool)]
    ax_f.scatter(
        all_meshes["lon"],
        all_meshes["lat"],
        s=4.5,
        color="#C9C9C9",
        alpha=0.42,
        linewidths=0,
        label="sampled meshes",
    )
    ax_f.scatter(
        tail["lon"],
        tail["lat"],
        s=6.0,
        color="#C44E52",
        alpha=0.50,
        linewidths=0,
        label="residual-tail meshes",
    )
    ax_f.scatter(
        ne_tail["lon"],
        ne_tail["lat"],
        s=8.5,
        color="#E58606",
        alpha=0.78,
        linewidths=0,
        label="NE-Hokkaido tails",
    )
    ne = all_meshes[all_meshes["is_ne_hokkaido"].astype(bool)]
    if not ne.empty:
        pad_lon = 0.45
        pad_lat = 0.30
        ax_f.add_patch(
            Rectangle(
                (ne["lon"].min() - pad_lon, ne["lat"].min() - pad_lat),
                ne["lon"].max() - ne["lon"].min() + 2 * pad_lon,
                ne["lat"].max() - ne["lat"].min() + 2 * pad_lat,
                fill=False,
                edgecolor="#E58606",
                linewidth=1.2,
            )
        )
    groups = physical_summary.set_index("group")
    ne_tail_pct = float(groups.loc["northeast_hokkaido_band", "mesh_tail_pct"])
    outside_tail_pct = float(groups.loc["outside_northeast_hokkaido", "mesh_tail_pct"])
    ne_d1400 = float(groups.loc["northeast_hokkaido_band", "nearest_d1400_median_m"])
    outside_d1400 = float(groups.loc["outside_northeast_hokkaido", "nearest_d1400_median_m"])
    note = (
        f"tail: {ne_tail_pct:.1f}% vs {outside_tail_pct:.1f}%\n"
        f"D1400: {ne_d1400:.0f} vs {outside_d1400:.0f} m"
    )
    ax_f.text(
        0.02,
        0.03,
        note,
        transform=ax_f.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#DDDDDD", alpha=0.90),
    )
    ax_f.set_xlim(128.0, 146.5)
    ax_f.set_ylim(30.0, 46.3)
    ax_f.set_xlabel("Longitude")
    ax_f.set_ylabel("Latitude")
    panel_title(ax_f, "F", "Northeast-Hokkaido residual-tail case")
    ax_f.grid(color="#EAEAEA", lw=0.7)
    ax_f.legend(frameon=False, loc="upper left", markerscale=1.8)

    fig.savefig(OUTDIR / "figure8.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure8.png", dpi=240, bbox_inches="tight")
    save_axes_subset(fig, [ax_a, ax_b, ax_c, ax_d], OUTDIR / "figure8a.pdf")
    save_axes_subset(fig, [ax_a, ax_b, ax_c, ax_d], OUTDIR / "figure8a.png")
    ax_e.texts[0].set_text("A")
    ax_f.texts[1].set_text("B")
    save_axes_subset(fig, [ax_e, ax_f], OUTDIR / "figure8b.pdf")
    save_axes_subset(fig, [ax_e, ax_f], OUTDIR / "figure8b.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
