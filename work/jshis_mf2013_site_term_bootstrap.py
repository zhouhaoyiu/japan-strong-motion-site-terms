#!/usr/bin/env python3
"""Cluster bootstrap uncertainty for MF2013 site-term improvements.

The official residual script writes event- and site-level aggregates. Those
aggregates are sufficient for clustered bootstrap confidence intervals of MAE
reductions without storing all per-record residuals.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVENT_TERMS = PROJECT_ROOT / "outputs" / "jshis_mf2013_official_event_residual_terms.csv"
DEFAULT_STATION_TERMS = PROJECT_ROOT / "outputs" / "jshis_mf2013_official_station_residual_terms.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"

TARGET_ORDER = ["PGA RotD50", "SA(0.3s) RotD50", "SA(1.0s) RotD50", "SA(3.0s) RotD50"]
NETWORK_LABELS = {1: "K-NET", 2: "KiK-net"}
SOURCE_CLASS_LABELS = {1: "CRUSTAL", 2: "INTERPLATE", 3: "INTRAPLATE"}


def paired_cluster_frame(
    frame: pd.DataFrame,
    cluster_col: str,
    filter_query: str | None = None,
) -> pd.DataFrame:
    work = frame.copy()
    if filter_query:
        work = work.query(filter_query).copy()
    keep = [
        "target",
        "target_label",
        "model",
        cluster_col,
        "n_records",
        "absolute_error_sum",
    ]
    for optional in ["network_label", "source_class_label"]:
        if optional in work.columns:
            keep.append(optional)
    work = work[keep].copy()
    rows = []
    for (target, label), sub in work.groupby(["target", "target_label"], dropna=False):
        basic = sub[sub["model"].eq("mf2013_basic")].set_index(cluster_col)
        site = sub[sub["model"].eq("mf2013_site")].set_index(cluster_col)
        common = basic.index.intersection(site.index)
        if common.empty:
            continue
        paired = pd.DataFrame(
            {
                "target": target,
                "target_label": label,
                cluster_col: common,
                "basic_n_records": basic.loc[common, "n_records"].to_numpy(dtype=float),
                "basic_absolute_error_sum": basic.loc[common, "absolute_error_sum"].to_numpy(dtype=float),
                "site_n_records": site.loc[common, "n_records"].to_numpy(dtype=float),
                "site_absolute_error_sum": site.loc[common, "absolute_error_sum"].to_numpy(dtype=float),
            }
        )
        rows.append(paired)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def bootstrap_one(
    paired: pd.DataFrame,
    cluster_col: str,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> dict[str, float]:
    basic_abs = paired["basic_absolute_error_sum"].to_numpy(dtype=float)
    basic_n = paired["basic_n_records"].to_numpy(dtype=float)
    site_abs = paired["site_absolute_error_sum"].to_numpy(dtype=float)
    site_n = paired["site_n_records"].to_numpy(dtype=float)

    basic_mae = float(basic_abs.sum() / basic_n.sum())
    site_mae = float(site_abs.sum() / site_n.sum())
    reduction = basic_mae - site_mae
    reduction_pct = reduction / basic_mae * 100.0

    n_clusters = len(paired)
    boot_abs = np.empty(n_bootstrap)
    boot_pct = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n_clusters, size=n_clusters)
        b_mae = basic_abs[idx].sum() / basic_n[idx].sum()
        s_mae = site_abs[idx].sum() / site_n[idx].sum()
        diff = b_mae - s_mae
        boot_abs[i] = diff
        boot_pct[i] = diff / b_mae * 100.0

    return {
        "n_clusters": int(n_clusters),
        "basic_records": int(basic_n.sum()),
        "site_records": int(site_n.sum()),
        "basic_mae": basic_mae,
        "site_mae": site_mae,
        "mae_reduction": reduction,
        "percent_mae_reduction": reduction_pct,
        "mae_reduction_ci_low": float(np.quantile(boot_abs, 0.025)),
        "mae_reduction_ci_high": float(np.quantile(boot_abs, 0.975)),
        "percent_reduction_ci_low": float(np.quantile(boot_pct, 0.025)),
        "percent_reduction_ci_high": float(np.quantile(boot_pct, 0.975)),
        "bootstrap_mean_percent_reduction": float(np.mean(boot_pct)),
        "p_reduction_le_zero": float(np.mean(boot_abs <= 0.0)),
    }


def bootstrap_table(
    frame: pd.DataFrame,
    cluster_col: str,
    cluster_level: str,
    n_bootstrap: int,
    rng: np.random.Generator,
    filter_query: str | None = None,
    subgroup: str = "all",
) -> pd.DataFrame:
    paired = paired_cluster_frame(frame, cluster_col, filter_query=filter_query)
    if paired.empty:
        return pd.DataFrame()
    rows = []
    for (target, target_label), sub in paired.groupby(["target", "target_label"]):
        item = bootstrap_one(sub, cluster_col, rng, n_bootstrap)
        item.update(
            {
                "cluster_level": cluster_level,
                "subgroup": subgroup,
                "target": target,
                "target_label": target_label,
                "cluster_column": cluster_col,
                "n_bootstrap": n_bootstrap,
            }
        )
        rows.append(item)
    return pd.DataFrame(rows)


def run(event_terms: Path, station_terms: Path, output_dir: Path, n_bootstrap: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    event = pd.read_csv(event_terms)
    station = pd.read_csv(station_terms)
    if "network_label" not in station.columns and "obs_network_id" in station.columns:
        station["network_label"] = station["obs_network_id"].map(NETWORK_LABELS)
    if "source_class_label" not in event.columns and "eq_location_type_id" in event.columns:
        event["source_class_label"] = event["eq_location_type_id"].map(SOURCE_CLASS_LABELS)
    rng = np.random.default_rng(20260605)

    pieces = [
        bootstrap_table(event, "eq_source_id", "event", n_bootstrap, rng),
        bootstrap_table(station, "siteid2", "site", n_bootstrap, rng),
    ]
    if "network_label" in station.columns:
        for network in sorted(station["network_label"].dropna().unique()):
            pieces.append(
                bootstrap_table(
                    station,
                    "siteid2",
                    "site",
                    n_bootstrap,
                    rng,
                    filter_query=f"network_label == '{network}'",
                    subgroup=str(network),
                )
            )
    if "source_class_label" in event.columns:
        for source_class in sorted(event["source_class_label"].dropna().unique()):
            pieces.append(
                bootstrap_table(
                    event,
                    "eq_source_id",
                    "event",
                    n_bootstrap,
                    rng,
                    filter_query=f"source_class_label == '{source_class}'",
                    subgroup=str(source_class),
                )
            )
    ci = pd.concat([part for part in pieces if not part.empty], ignore_index=True)
    ci["target_order"] = ci["target_label"].map({label: idx for idx, label in enumerate(TARGET_ORDER)})
    ci = ci.sort_values(["cluster_level", "subgroup", "target_order"]).drop(columns=["target_order"])
    ci.to_csv(output_dir / "jshis_mf2013_site_term_bootstrap_ci.csv", index=False)
    plot_ci(ci, output_dir)
    write_report(ci, output_dir)


def plot_ci(ci: pd.DataFrame, output_dir: Path) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    sub = ci[ci["subgroup"].eq("all")].copy()
    if sub.empty:
        return
    sub["target_label"] = pd.Categorical(sub["target_label"], categories=TARGET_ORDER, ordered=True)
    sub = sub.sort_values(["cluster_level", "target_label"])
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    offsets = {"event": -0.11, "site": 0.11}
    colors = {"event": "#1f77b4", "site": "#ff7f0e"}
    x_base = np.arange(len(TARGET_ORDER))
    for level, group in sub.groupby("cluster_level"):
        x = np.array([TARGET_ORDER.index(label) for label in group["target_label"]]) + offsets.get(level, 0.0)
        y = group["percent_mae_reduction"].to_numpy()
        yerr = np.vstack(
            [
                y - group["percent_reduction_ci_low"].to_numpy(),
                group["percent_reduction_ci_high"].to_numpy() - y,
            ]
        )
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            capsize=4,
            label=f"{level}-cluster bootstrap",
            color=colors.get(level),
        )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x_base)
    ax.set_xticklabels(TARGET_ORDER, rotation=15, ha="right")
    ax.set_ylabel("MF2013 site-term MAE reduction (%)")
    ax.set_title("Cluster Bootstrap CIs for Official MF2013 Site-Term Gains")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "jshis_mf2013_site_term_bootstrap_ci.png", dpi=180)
    plt.close(fig)


def write_report(ci: pd.DataFrame, output_dir: Path) -> None:
    lines = ["# MF2013 Site-Term Bootstrap Uncertainty\n\n"]
    lines.append("## Method\n")
    lines.append("- Inputs: official J-SHIS MF2013 event and station residual aggregate tables.\n")
    lines.append("- Comparison: `mf2013_site` minus `mf2013_basic`, using paired clusters present in both models.\n")
    lines.append("- Statistic: MAE reduction and percent MAE reduction from official `D1400`/`AVS30` site terms.\n")
    n_boot = int(ci["n_bootstrap"].iloc[0]) if not ci.empty else 0
    lines.append(f"- Bootstrap replicates per row: {n_boot:,}.\n\n")

    main = ci[ci["subgroup"].eq("all")]
    lines.append("## Overall clustered intervals\n")
    for cluster_level in ["event", "site"]:
        lines.append(f"### {cluster_level}-cluster bootstrap\n")
        sub = main[main["cluster_level"].eq(cluster_level)].copy()
        sub["target_label"] = pd.Categorical(sub["target_label"], categories=TARGET_ORDER, ordered=True)
        for _, row in sub.sort_values("target_label").iterrows():
            lines.append(
                f"- {row['target_label']}: reduction={row['percent_mae_reduction']:.1f}% "
                f"(95% CI {row['percent_reduction_ci_low']:.1f} to {row['percent_reduction_ci_high']:.1f}%), "
                f"p(reduction<=0)={row['p_reduction_le_zero']:.4f}, clusters={int(row['n_clusters'])}.\n"
            )
        lines.append("\n")

    lines.append("## Network and source-class sensitivity\n")
    sensitivity = ci[~ci["subgroup"].eq("all")].copy()
    if sensitivity.empty:
        lines.append("- No subgroup sensitivity rows were generated.\n")
    else:
        for subgroup in sorted(sensitivity["subgroup"].unique()):
            sub = sensitivity[sensitivity["subgroup"].eq(subgroup)].copy()
            best = sub[sub["target_label"].isin(["SA(1.0s) RotD50", "SA(3.0s) RotD50"])]
            if best.empty:
                continue
            vals = "; ".join(
                f"{row['target_label']} {row['percent_mae_reduction']:.1f}% "
                f"[{row['percent_reduction_ci_low']:.1f}, {row['percent_reduction_ci_high']:.1f}]"
                for _, row in best.iterrows()
            )
            lines.append(f"- {subgroup}: {vals}.\n")

    lines.append("\n## Outputs\n")
    lines.append("- `jshis_mf2013_site_term_bootstrap_ci.csv`\n")
    lines.append("- `figures/jshis_mf2013_site_term_bootstrap_ci.png`\n")
    lines.append("\n## Guardrail\n")
    lines.append(
        "The intervals are clustered over events or sites using aggregate absolute-error sums. "
        "They quantify robustness of MAE reductions, not uncertainty in the MF2013 coefficients themselves.\n"
    )
    (output_dir / "jshis_mf2013_site_term_bootstrap.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--event-terms", type=Path, default=DEFAULT_EVENT_TERMS)
    parser.add_argument("--station-terms", type=Path, default=DEFAULT_STATION_TERMS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    args = parser.parse_args()
    run(args.event_terms, args.station_terms, args.output_dir, args.n_bootstrap)


if __name__ == "__main__":
    main()
