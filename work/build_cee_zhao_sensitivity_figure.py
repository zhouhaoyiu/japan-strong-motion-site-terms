#!/usr/bin/env python3
"""Build a CEE-style external-GMPE sensitivity figure from Zhao 2006 audit outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "outputs"
FIGURE_DIR = OUTPUT_DIR / "figures"

TARGET_ORDER = ["PGA RotD50", "SA(0.3s) RotD50", "SA(1.0s) RotD50", "SA(3.0s) RotD50"]
TARGET_SHORT = {
    "PGA RotD50": "PGA",
    "SA(0.3s) RotD50": "SA 0.3 s",
    "SA(1.0s) RotD50": "SA 1.0 s",
    "SA(3.0s) RotD50": "SA 3.0 s",
}
FEATURE_ORDER = ["log10(VS20)", "log10(AVS30)", "log10(D1100)", "log10(D1400)", "log10(Dbase)"]
FEATURE_SHORT = {
    "log10(VS20)": "VS20",
    "log10(AVS30)": "AVS30",
    "log10(D1100)": "D1100",
    "log10(D1400)": "D1400",
    "log10(Dbase)": "Dbase",
}
SOURCE_ORDER = ["CRUSTAL", "INTERPLATE", "INTRAPLATE"]
SOURCE_COLORS = {
    "CRUSTAL": "#0072B2",
    "INTERPLATE": "#D55E00",
    "INTRAPLATE": "#009E73",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{stem}.png")
    fig.savefig(FIGURE_DIR / f"{stem}.pdf")


def d1400_model_sensitivity(ax: plt.Axes, mf: pd.DataFrame, zhao: pd.DataFrame) -> pd.DataFrame:
    mf_d = mf[(mf["model"].eq("mf2013_basic")) & (mf["feature"].eq("log10(D1400)"))].copy()
    zhao_d = zhao[zhao["feature"].eq("log10(D1400)")].copy()
    mf_d["model_family"] = "MF2013 basic"
    zhao_d["model_family"] = "Zhao 2006 audit"
    combined = pd.concat(
        [
            mf_d[["target_label", "model_family", "spearman_rho", "n_sites"]],
            zhao_d[["target_label", "model_family", "spearman_rho", "n_sites"]],
        ],
        ignore_index=True,
    )
    combined["target_label"] = pd.Categorical(combined["target_label"], categories=TARGET_ORDER, ordered=True)
    combined = combined.sort_values(["model_family", "target_label"])

    x = np.arange(len(TARGET_ORDER))
    colors = {"MF2013 basic": "#2166ac", "Zhao 2006 audit": "#b2182b"}
    markers = {"MF2013 basic": "o", "Zhao 2006 audit": "s"}
    for name, sub in combined.groupby("model_family", sort=False):
        y = [sub[sub["target_label"].eq(label)]["spearman_rho"].iloc[0] for label in TARGET_ORDER]
        ax.plot(x, y, marker=markers[name], ms=5.5, lw=1.8, label=name, color=colors[name])
    ax.axhline(0, color="#333333", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([TARGET_SHORT[label] for label in TARGET_ORDER])
    ax.set_ylim(-0.10, 0.70)
    ax.set_ylabel("Spearman rho with log10(D1400)")
    ax.set_title("A  Long-period basin signal is preserved")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False, loc="upper left")
    return combined


def zhao_heatmap(ax: plt.Axes, zhao: pd.DataFrame) -> pd.DataFrame:
    sub = zhao[zhao["feature"].isin(FEATURE_ORDER)].copy()
    sub["target_label"] = pd.Categorical(sub["target_label"], categories=TARGET_ORDER, ordered=True)
    sub["feature"] = pd.Categorical(sub["feature"], categories=FEATURE_ORDER, ordered=True)
    heat = sub.pivot(index="target_label", columns="feature", values="spearman_rho").loc[TARGET_ORDER, FEATURE_ORDER]
    im = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-0.55, vmax=0.55, aspect="auto")
    ax.set_yticks(np.arange(len(TARGET_ORDER)))
    ax.set_yticklabels([TARGET_SHORT[label] for label in TARGET_ORDER])
    ax.set_xticks(np.arange(len(FEATURE_ORDER)))
    ax.set_xticklabels([FEATURE_SHORT[label] for label in FEATURE_ORDER], rotation=35, ha="right")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.5, color="white" if abs(val) >= 0.35 else "#111111")
    ax.set_title("B  Zhao residual-site associations")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.035)
    cbar.set_label("rho", labelpad=6)
    return heat.reset_index()


def source_coverage(ax: plt.Axes, by_source: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    source = by_source[by_source["target_label"].eq("PGA RotD50")].copy()
    source = source.set_index("source_class_label").loc[SOURCE_ORDER].reset_index()
    x = np.arange(len(source))
    counts = source["n_records"].to_numpy()
    bars = ax.bar(x, counts / 1000.0, color=[SOURCE_COLORS[label] for label in source["source_class_label"]])
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2.2, f"{int(count):,}", ha="center", va="bottom", fontsize=8)
    usable = int(coverage["usable_rows"].iloc[0])
    invalid = int(coverage["invalid_rows"].iloc[0])
    ax.set_xticks(x)
    ax.set_xticklabels(source["source_class_label"], rotation=20, ha="right")
    ax.set_ylabel("Records (k)", labelpad=8)
    ax.set_title("C  External-GMPE audit coverage")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.text(
        0.02,
        0.08,
        f"{usable:,} usable rows per target\n{invalid:,} missing/invalid rows",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc"),
    )
    return source[["source_class_label", "n_records"]]


def write_source_table(
    d1400: pd.DataFrame,
    heat: pd.DataFrame,
    source: pd.DataFrame,
    output_path: Path,
) -> None:
    rows = []
    for _, row in d1400.iterrows():
        rows.append(
            {
                "panel": "A",
                "quantity": "d1400_model_sensitivity",
                "target_label": row["target_label"],
                "feature": "log10(D1400)",
                "model_family": row["model_family"],
                "value": row["spearman_rho"],
                "n": row["n_sites"],
            }
        )
    for _, row in heat.iterrows():
        target = row["target_label"]
        for feature in FEATURE_ORDER:
            rows.append(
                {
                    "panel": "B",
                    "quantity": "zhao_site_correlation",
                    "target_label": target,
                    "feature": feature,
                    "model_family": "Zhao 2006 audit",
                    "value": row[feature],
                    "n": np.nan,
                }
            )
    for _, row in source.iterrows():
        rows.append(
            {
                "panel": "C",
                "quantity": "source_class_coverage",
                "target_label": "all targets",
                "feature": row["source_class_label"],
                "model_family": "Zhao 2006 audit",
                "value": row["n_records"],
                "n": row["n_records"],
            }
        )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def main() -> None:
    setup_style()
    mf = pd.read_csv(OUTPUT_DIR / "jshis_mf2013_official_station_site_correlations.csv")
    zhao = pd.read_csv(OUTPUT_DIR / "jshis_zhao2006_external_gmpe_station_site_correlations.csv")
    by_source = pd.read_csv(OUTPUT_DIR / "jshis_zhao2006_external_gmpe_by_source_class_metrics.csv")
    coverage = pd.read_csv(OUTPUT_DIR / "jshis_zhao2006_external_gmpe_coverage_summary.csv")

    fig = plt.figure(figsize=(13.2, 5.4))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.42, 1.00], wspace=0.52)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    d1400 = d1400_model_sensitivity(ax_a, mf, zhao)
    heat = zhao_heatmap(ax_b, zhao)
    source = source_coverage(ax_c, by_source, coverage)

    fig.text(
        0.01,
        0.01,
        "Zhao 2006 audit uses OpenQuake hazardlib, J-SHIS fault_dist as rrup proxy, and vs30/avs30 as a site proxy.",
        fontsize=8,
        color="#555555",
    )
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.17, top=0.88)
    save_figure(fig, "cee_fig5_zhao2006_external_gmpe_sensitivity")
    plt.close(fig)

    write_source_table(
        d1400=d1400,
        heat=heat,
        source=source,
        output_path=OUTPUT_DIR / "cee_fig5_zhao2006_external_gmpe_sensitivity_source.csv",
    )


if __name__ == "__main__":
    main()
