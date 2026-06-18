#!/usr/bin/env python3
"""Quantify SA(3.0 s) station-correction impacts on hazard curves and UHS.

This script uses out-of-fold spatial-block station-residual predictions from
the station model. It does not run a full PSHA calculation because source
recurrence rates and official hazard curves are not part of the local inputs.
Instead, it applies the exact curve-shift relationship implied by a station
correction in log10 spectral acceleration.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PREDICTIONS = PROJECT_ROOT / "outputs" / "jshis_station_model_ablation_predictions.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures"

SA3_LABEL = "SA(3.0s) RotD50"


def weighted_quantile(values: np.ndarray, quantiles: list[float], weights: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    quantiles_arr = np.asarray(quantiles, dtype=float)
    if weights is None:
        return np.quantile(values, quantiles_arr)
    weights = np.asarray(weights, dtype=float)
    sorter = np.argsort(values)
    values = values[sorter]
    weights = weights[sorter]
    cdf = np.cumsum(weights) - 0.5 * weights
    cdf = cdf / np.sum(weights)
    return np.interp(quantiles_arr, cdf, values)


def nearest_quantile_rows(frame: pd.DataFrame, quantiles: list[float]) -> pd.DataFrame:
    rows = []
    factor = frame["sa3_uhs_multiplier"].to_numpy(dtype=float)
    for q in quantiles:
        value = float(np.quantile(factor, q))
        idx = int(np.argmin(np.abs(factor - value)))
        row = frame.iloc[idx].copy()
        row["example_type"] = f"q{int(q * 100):02d}"
        row["target_quantile_factor"] = value
        rows.append(row)
    high = frame.nlargest(1, "sa3_uhs_multiplier").iloc[0].copy()
    high["example_type"] = "maximum_positive"
    high["target_quantile_factor"] = high["sa3_uhs_multiplier"]
    low = frame.nsmallest(1, "sa3_uhs_multiplier").iloc[0].copy()
    low["example_type"] = "maximum_negative"
    low["target_quantile_factor"] = low["sa3_uhs_multiplier"]
    rows.extend([high, low])
    out = pd.DataFrame(rows)
    preferred = [
        "example_type",
        "site_code",
        "network_label",
        "n_records",
        "mean_residual",
        "sa3_correction_log10",
        "sa3_uhs_multiplier",
        "fixed_sa_hazard_ratio_beta2",
        "fixed_sa_hazard_ratio_beta3",
        "fixed_sa_hazard_ratio_beta4",
        "corrected_sa3_if_baseline_0p20g",
        "corrected_sa3_if_baseline_0p35g",
    ]
    return out[preferred].reset_index(drop=True)


def build_station_impact(predictions_path: Path) -> pd.DataFrame:
    pred = pd.read_csv(predictions_path)
    sa3 = pred[pred["split"].eq("spatial_block") & pred["target_label"].astype(str).eq(SA3_LABEL)].copy()
    if sa3.empty:
        raise ValueError(f"No spatial-block predictions found for {SA3_LABEL}")
    sa3 = sa3.rename(columns={"full_site_location_prediction": "sa3_correction_log10"})
    sa3["sa3_uhs_multiplier"] = np.power(10.0, sa3["sa3_correction_log10"].astype(float))
    sa3["sa3_uhs_delta_pct"] = 100.0 * (sa3["sa3_uhs_multiplier"] - 1.0)
    for beta in [2, 3, 4]:
        sa3[f"fixed_sa_hazard_ratio_beta{beta}"] = np.power(sa3["sa3_uhs_multiplier"], beta)
        sa3[f"fixed_sa_hazard_delta_pct_beta{beta}"] = 100.0 * (sa3[f"fixed_sa_hazard_ratio_beta{beta}"] - 1.0)
    # Absolute ordinates are illustrative scalings. The multipliers are the
    # actual result and can be applied to any baseline UHS ordinate.
    sa3["corrected_sa3_if_baseline_0p20g"] = 0.20 * sa3["sa3_uhs_multiplier"]
    sa3["corrected_sa3_if_baseline_0p35g"] = 0.35 * sa3["sa3_uhs_multiplier"]
    keep = [
        "siteid2",
        "site_code",
        "network_label",
        "n_records",
        "mean_residual",
        "fold",
        "sa3_correction_log10",
        "sa3_uhs_multiplier",
        "sa3_uhs_delta_pct",
        "fixed_sa_hazard_ratio_beta2",
        "fixed_sa_hazard_delta_pct_beta2",
        "fixed_sa_hazard_ratio_beta3",
        "fixed_sa_hazard_delta_pct_beta3",
        "fixed_sa_hazard_ratio_beta4",
        "fixed_sa_hazard_delta_pct_beta4",
        "corrected_sa3_if_baseline_0p20g",
        "corrected_sa3_if_baseline_0p35g",
    ]
    return sa3[keep].sort_values("sa3_uhs_multiplier").reset_index(drop=True)


def build_summary(impact: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "sa3_correction_log10",
        "sa3_uhs_multiplier",
        "sa3_uhs_delta_pct",
        "fixed_sa_hazard_ratio_beta2",
        "fixed_sa_hazard_ratio_beta3",
        "fixed_sa_hazard_ratio_beta4",
    ]
    quantiles = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
    rows = []
    weights = impact["n_records"].to_numpy(dtype=float)
    for metric in metrics:
        values = impact[metric].to_numpy(dtype=float)
        row = {
            "metric": metric,
            "n_stations": len(impact),
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
        }
        for q, value in zip(quantiles, np.quantile(values, quantiles), strict=True):
            row[f"q{int(q * 100):02d}"] = float(value)
        wq = weighted_quantile(values, [0.05, 0.50, 0.95], weights)
        row["record_weighted_q05"] = float(wq[0])
        row["record_weighted_q50"] = float(wq[1])
        row["record_weighted_q95"] = float(wq[2])
        rows.append(row)
    return pd.DataFrame(rows)


def build_hazard_curve_points(examples: pd.DataFrame, beta: float = 3.0) -> pd.DataFrame:
    x = np.geomspace(0.15, 3.0, 180)
    rows = []
    for _, row in examples.iterrows():
        factor = float(row["sa3_uhs_multiplier"])
        for xx in x:
            rows.append(
                {
                    "example_type": row["example_type"],
                    "site_code": row["site_code"],
                    "network_label": row["network_label"],
                    "hazard_slope_beta": beta,
                    "sa3_over_baseline_475yr": xx,
                    "baseline_exceedance_ratio": xx ** (-beta),
                    "corrected_exceedance_ratio": (xx / factor) ** (-beta),
                    "sa3_uhs_multiplier": factor,
                }
            )
    return pd.DataFrame(rows)


def write_report(summary: pd.DataFrame, examples: pd.DataFrame, output_path: Path) -> None:
    factor = summary[summary["metric"].eq("sa3_uhs_multiplier")].iloc[0]
    hazard = summary[summary["metric"].eq("fixed_sa_hazard_ratio_beta3")].iloc[0]
    lines = [
        "# SA(3.0 s) Hazard-Curve and UHS Impact",
        "",
        "This analysis is a hazard-impact transformation, not a full PSHA run. The local inputs do not include source recurrence rates or official site-specific hazard curves.",
        "",
        "If the station correction is c in log10 spectral acceleration, the multiplicative SA(3.0 s) factor is f = 10^c. For any baseline hazard curve lambda0(S), the corrected curve is lambda_c(S) = lambda0(S / f). At a fixed annual exceedance rate, the SA(3.0 s) UHS ordinate becomes f times the baseline ordinate. For a local power-law hazard slope lambda0(S) proportional to S^-beta, the fixed-ordinate annual exceedance-rate ratio is f^beta.",
        "",
        "## Main station-level distribution",
        "",
        f"- Stations: {int(factor['n_stations'])}",
        f"- SA(3.0 s) UHS multiplier q05/q50/q95: {factor['q05']:.3f} / {factor['q50']:.3f} / {factor['q95']:.3f}",
        f"- SA(3.0 s) UHS multiplier range: {factor['min']:.3f} to {factor['max']:.3f}",
        f"- Fixed-SA hazard-rate ratio with beta=3 q05/q50/q95: {hazard['q05']:.3f} / {hazard['q50']:.3f} / {hazard['q95']:.3f}",
        "",
        "## Example stations",
        "",
        "| example_type | site_code | network_label | correction_log10 | UHS_multiplier | hazard_ratio_beta3 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in examples.iterrows():
        lines.append(
            f"| {row['example_type']} | {row['site_code']} | {row['network_label']} | "
            f"{row['sa3_correction_log10']:.3f} | {row['sa3_uhs_multiplier']:.3f} | "
            f"{row['fixed_sa_hazard_ratio_beta3']:.3f} |"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_figure(impact: pd.DataFrame, examples: pd.DataFrame, curves: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig = plt.figure(figsize=(12.6, 8.6))
    grid = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    factors = impact["sa3_uhs_multiplier"].to_numpy(dtype=float)
    ax_a.hist(factors, bins=np.geomspace(max(factors.min(), 0.08), factors.max() * 1.05, 34), color="#4C78A8", edgecolor="white")
    ax_a.set_xscale("log")
    ax_a.axvline(1.0, color="black", lw=1.0)
    for q, ls in [(0.05, "--"), (0.50, "-"), (0.95, "--")]:
        ax_a.axvline(np.quantile(factors, q), color="#E45756", lw=1.0, ls=ls)
    ax_a.set_xlabel("SA(3.0 s) UHS multiplier")
    ax_a.set_ylabel("Number of stations")
    ax_a.set_title("A Distribution of SA(3.0 s) station-correction factors")
    ax_a.grid(axis="y", color="#e6e6e6")

    color_map = {
        "maximum_negative": "#4C78A8",
        "q05": "#72B7B2",
        "q50": "#777777",
        "q95": "#F58518",
        "maximum_positive": "#E45756",
    }
    selected = examples[examples["example_type"].isin(["maximum_negative", "q50", "maximum_positive"])]
    x = curves["sa3_over_baseline_475yr"].drop_duplicates().to_numpy(dtype=float)
    ax_b.plot(x, x ** -3.0, color="black", lw=1.4, label="baseline")
    for _, row in selected.iterrows():
        sub = curves[curves["example_type"].eq(row["example_type"])]
        ax_b.plot(
            sub["sa3_over_baseline_475yr"],
            sub["corrected_exceedance_ratio"],
            lw=1.4,
            color=color_map.get(row["example_type"], "#999999"),
            label=f"{row['example_type']} ({row['site_code']}, x{row['sa3_uhs_multiplier']:.2f})",
        )
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlabel("SA(3.0 s) / baseline 475-year SA")
    ax_b.set_ylabel("Annual exceedance ratio")
    ax_b.set_title("B Hazard-curve shift for local slope beta=3")
    ax_b.legend(frameon=False, fontsize=7.5)
    ax_b.grid(color="#e6e6e6", which="both")

    q_examples = examples[examples["example_type"].isin(["q05", "q25", "q50", "q75", "q95"])].copy()
    ax_c.bar(
        np.arange(len(q_examples)),
        q_examples["sa3_uhs_multiplier"],
        color=["#4C78A8", "#72B7B2", "#777777", "#F58518", "#E45756"],
        edgecolor="black",
        linewidth=0.5,
    )
    ax_c.axhline(1.0, color="black", lw=1.0)
    ax_c.set_xticks(np.arange(len(q_examples)))
    ax_c.set_xticklabels(q_examples["example_type"])
    ax_c.set_ylabel("Corrected / baseline SA(3.0 s) UHS")
    ax_c.set_title("C UHS ordinate ratio at any return period")
    ax_c.grid(axis="y", color="#e6e6e6")

    slopes = np.array([2, 3, 4], dtype=float)
    for _, row in q_examples.iterrows():
        y = np.power(row["sa3_uhs_multiplier"], slopes)
        ax_d.plot(slopes, y, marker="o", lw=1.4, label=f"{row['example_type']} x{row['sa3_uhs_multiplier']:.2f}")
    ax_d.axhline(1.0, color="black", lw=1.0)
    ax_d.set_xlabel("Local hazard slope beta")
    ax_d.set_ylabel("Fixed-SA exceedance-rate ratio")
    ax_d.set_title("D Sensitivity of hazard-rate ratio to curve slope")
    ax_d.set_xticks(slopes)
    ax_d.set_yscale("log")
    ax_d.legend(frameon=False, fontsize=7.5, ncols=2)
    ax_d.grid(color="#e6e6e6", which="both")

    fig.suptitle("Effect of SA(3.0 s) station corrections on hazard curves and UHS ordinates", fontsize=13)
    fig.savefig(figure_dir / "cee_fig8_hazard_impact.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / "cee_fig8_hazard_impact.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    impact = build_station_impact(args.predictions)
    summary = build_summary(impact)
    examples = nearest_quantile_rows(impact, [0.05, 0.25, 0.50, 0.75, 0.95])
    curves = build_hazard_curve_points(examples)

    impact.to_csv(args.output_dir / "jshis_station_hazard_impact_by_station.csv", index=False)
    summary.to_csv(args.output_dir / "jshis_station_hazard_impact_summary.csv", index=False)
    examples.to_csv(args.output_dir / "jshis_station_hazard_impact_examples.csv", index=False)
    curves.to_csv(args.output_dir / "jshis_station_hazard_curve_points.csv", index=False)
    write_report(summary, examples, args.output_dir / "jshis_station_hazard_impact.md")
    build_figure(impact, examples, curves, args.figure_dir)

    print(f"Wrote {args.output_dir / 'jshis_station_hazard_impact_by_station.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_hazard_impact_summary.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_hazard_impact_examples.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_hazard_curve_points.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_station_hazard_impact.md'}")
    print(f"Wrote {args.figure_dir / 'cee_fig8_hazard_impact.pdf'}")


if __name__ == "__main__":
    main()
