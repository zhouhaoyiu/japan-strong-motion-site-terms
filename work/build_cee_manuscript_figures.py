#!/usr/bin/env python3
"""Build manuscript-style figure drafts for the CEE submission direction."""

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


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURE_DIR / f"{stem}.png")
    fig.savefig(FIGURE_DIR / f"{stem}.pdf")


def fig1_workflow() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    ax.axis("off")

    boxes = [
        (0.05, 0.66, 0.23, 0.20, "Official J-SHIS/NIED\nflatfile sub1-v2024\n333,808 records"),
        (0.38, 0.66, 0.23, 0.20, "MF2013 core equation\nMw, fault distance,\nsource class"),
        (0.72, 0.66, 0.23, 0.20, "D1400 / AVS30\nsite correction\nofficial site schema"),
        (0.05, 0.28, 0.23, 0.20, "Residual partition\nby target, event,\nsite, network"),
        (0.38, 0.28, 0.23, 0.20, "Cluster bootstrap\n1,737 events\n2,571 sites"),
        (0.72, 0.28, 0.23, 0.20, "Frequency-dependent\nsite-response evidence\nSA 1.0-3.0 s"),
    ]
    colors = ["#e8f1fb", "#f0f4e8", "#fff2df", "#f3edf7", "#e9f7f0", "#fdecea"]
    for (x, y, w, h, text), color in zip(boxes, colors):
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor="#4a4a4a", linewidth=1.1)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", linespacing=1.25)

    arrows = [
        ((0.28, 0.76), (0.38, 0.76)),
        ((0.61, 0.76), (0.72, 0.76)),
        ((0.84, 0.66), (0.84, 0.48)),
        ((0.72, 0.38), (0.61, 0.38)),
        ((0.38, 0.38), (0.28, 0.38)),
    ]
    for start, end in arrows:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(arrowstyle="->", lw=1.4, color="#333333"),
        )

    ax.text(0.05, 0.93, "A", fontsize=14, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.08,
        0.93,
        "Analysis workflow for official-data MF2013 residual evidence",
        fontsize=13,
        fontweight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0.05,
        0.08,
        "The paper backbone is official residual evidence; machine-learning source-path-site modeling is an extension layer.",
        fontsize=9,
        color="#555555",
        transform=ax.transAxes,
    )
    fig.tight_layout()
    save_figure(fig, "cee_fig1_workflow_data_evidence")
    plt.close(fig)


def fig2_gain_ci(ci: pd.DataFrame) -> None:
    main = ci[ci["subgroup"].eq("all")].copy()
    main["target_label"] = pd.Categorical(main["target_label"], categories=TARGET_ORDER, ordered=True)
    main = main.sort_values(["cluster_level", "target_label"])
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    offsets = {"event": -0.10, "site": 0.10}
    colors = {"event": "#2166ac", "site": "#b2182b"}
    x_base = np.arange(len(TARGET_ORDER))
    for level, sub in main.groupby("cluster_level"):
        x = np.array([TARGET_ORDER.index(label) for label in sub["target_label"]]) + offsets[level]
        y = sub["percent_mae_reduction"].to_numpy()
        yerr = np.vstack(
            [
                y - sub["percent_reduction_ci_low"].to_numpy(),
                sub["percent_reduction_ci_high"].to_numpy() - y,
            ]
        )
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            capsize=4,
            ms=6,
            lw=1.4,
            label=f"{level}-cluster CI",
            color=colors[level],
        )
    ax.axhline(0, color="#222222", lw=0.8)
    ax.set_xticks(x_base)
    ax.set_xticklabels([TARGET_SHORT[label] for label in TARGET_ORDER])
    ax.set_ylabel("MAE reduction from D1400/AVS30 site terms (%)")
    ax.set_title("MF2013 site terms improve response spectra more than PGA")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.text(-0.42, 23.0, "B", fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "cee_fig2_mf2013_site_term_gain_ci")
    plt.close(fig)


def fig3_correlation_shift(effects: pd.DataFrame) -> None:
    corr = effects[effects["effect_type"].eq("station_residual_spearman_absrho")].copy()
    selected = [
        ("SA(3.0s) RotD50", "log10(D1400)"),
        ("SA(3.0s) RotD50", "log10(Dbase)"),
        ("SA(3.0s) RotD50", "log10(AVS30)"),
        ("SA(1.0s) RotD50", "log10(VS20)"),
        ("SA(1.0s) RotD50", "log10(D1400)"),
        ("SA(0.3s) RotD50", "log10(VS20)"),
        ("PGA RotD50", "log10(VS20)"),
    ]
    rows = []
    for target, feature in selected:
        sub = corr[(corr["target_label"].eq(target)) & (corr["feature"].eq(feature))]
        if not sub.empty:
            rows.append(sub.iloc[0])
    plot = pd.DataFrame(rows)
    labels = [f"{TARGET_SHORT[row['target_label']]}\n{row['feature'].replace('log10(', '').replace(')', '')}" for _, row in plot.iterrows()]
    x = np.arange(len(plot))
    fig, ax = plt.subplots(figsize=(10.0, 5.2))
    for xi, row in zip(x, plot.itertuples()):
        ax.plot([xi, xi], [row.basic_value, row.site_value], color="#9e9e9e", lw=1.4, zorder=1)
    ax.scatter(x - 0.035, plot["basic_value"], color="#2166ac", s=48, label="MF2013 basic residual", zorder=2)
    ax.scatter(x + 0.035, plot["site_value"], color="#b2182b", s=48, label="After D1400/AVS30 correction", zorder=3)
    ax.axhline(0, color="#222222", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Spearman rho with station mean residual")
    ax.set_title("Official site terms remove much of the long-period basin residual structure")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.text(-0.55, 0.66, "C", fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "cee_fig3_residual_site_correlation_shift")
    plt.close(fig)


def fig4_subgroup_robustness(ci: pd.DataFrame) -> None:
    sub = ci[
        ci["subgroup"].isin(["K-NET", "KiK-net", "CRUSTAL", "INTERPLATE", "INTRAPLATE"])
        & ci["target_label"].isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])
    ].copy()
    subgroup_order = ["K-NET", "KiK-net", "CRUSTAL", "INTERPLATE", "INTRAPLATE"]
    target_order = ["SA(1.0s) RotD50", "SA(3.0s) RotD50"]
    pivot = sub.pivot(index="subgroup", columns="target_label", values="percent_mae_reduction")
    low = sub.pivot(index="subgroup", columns="target_label", values="percent_reduction_ci_low")
    high = sub.pivot(index="subgroup", columns="target_label", values="percent_reduction_ci_high")
    pivot = pivot.loc[subgroup_order, target_order]
    low = low.loc[subgroup_order, target_order]
    high = high.loc[subgroup_order, target_order]

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    im = ax.imshow(pivot.to_numpy(), cmap="YlGnBu", vmin=0, vmax=28, aspect="auto")
    ax.set_xticks(np.arange(len(target_order)))
    ax.set_xticklabels(["SA 1.0 s", "SA 3.0 s"])
    ax.set_yticks(np.arange(len(subgroup_order)))
    ax.set_yticklabels(subgroup_order)
    for i, group in enumerate(subgroup_order):
        for j, target in enumerate(target_order):
            val = pivot.loc[group, target]
            ax.text(
                j,
                i,
                f"{val:.1f}%\n[{low.loc[group, target]:.1f}, {high.loc[group, target]:.1f}]",
                ha="center",
                va="center",
                fontsize=8,
                color="white" if val >= 18 else "#111111",
            )
    ax.set_title("Response-spectrum gains are stable across networks and source classes")
    fig.colorbar(im, ax=ax, label="MAE reduction (%)", fraction=0.046, pad=0.04)
    ax.text(-0.75, -0.78, "D", fontsize=14, fontweight="bold", transform=ax.transData)
    fig.tight_layout()
    save_figure(fig, "cee_fig4_subgroup_robustness")
    plt.close(fig)


def write_caption_drafts() -> None:
    text = """# CEE Manuscript Figure Caption Drafts

## Figure 1

Workflow for the official-data residual analysis. The analysis links J-SHIS/NIED strong-motion records to official source and site metadata, evaluates MF2013 core and site-corrected predictions, partitions residuals by event and site, and tests whether the site-term gains are frequency dependent and robust under clustered resampling.

## Figure 2

MAE reduction from adding official MF2013 D1400/AVS30 site terms. Error bars are 95% cluster-bootstrap intervals over events and sites. Site terms provide weak improvement for PGA but large and stable improvement for 1.0 s and 3.0 s RotD50 spectral acceleration.

## Figure 3

Station residual correlations with official site parameters before and after MF2013 site correction. Basin-depth correlations at long periods and shallow-velocity correlations at short-to-intermediate periods are strong in the basic residuals; D1400/AVS30 corrections remove much of the long-period basin residual structure.

## Figure 4

Subgroup sensitivity of MF2013 site-term gains for SA(1.0s) and SA(3.0s). Values are percent MAE reductions with 95% cluster-bootstrap intervals. Positive gains across K-NET, KiK-net, and the three official source classes support a frequency-dependent site-response interpretation rather than a single-network artifact.

## Guardrail

These are figure drafts for manuscript planning. Final captions should include exact data version, residual definition, and limitations on omitted MF2013 AI/PH terms.
"""
    (OUTPUT_DIR / "cee_figure_caption_drafts.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_style()
    ci = pd.read_csv(OUTPUT_DIR / "jshis_mf2013_site_term_bootstrap_ci.csv")
    effects = pd.read_csv(OUTPUT_DIR / "jshis_mf2013_official_site_correction_effects.csv")
    fig1_workflow()
    fig2_gain_ci(ci)
    fig3_correlation_shift(effects)
    fig4_subgroup_robustness(ci)
    write_caption_drafts()


if __name__ == "__main__":
    main()
