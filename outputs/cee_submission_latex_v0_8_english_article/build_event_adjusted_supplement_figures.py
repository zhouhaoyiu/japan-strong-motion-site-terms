#!/usr/bin/env python3
"""Build supplementary figures from the event-adjusted release tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ARTICLE_DIR = Path(__file__).resolve().parent
SUPPLEMENT_DIR = ARTICLE_DIR / "supplement"
FIGURE_DIR = ARTICLE_DIR / "figures"

PERIOD_TICKS = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0]


def style() -> None:
    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def period_axis(ax: plt.Axes) -> None:
    ax.set_xscale("log")
    ax.set_xticks(PERIOD_TICKS, labels=["0.1", "0.2", "0.3", "0.5", "1", "2", "3", "5"])
    ax.set_xlabel("Period (s)")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURE_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


def spatial_validation_figure() -> None:
    predictions = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv")
    metrics = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_metrics.csv")
    frame = predictions[
        predictions["period_s"].eq(3.0) & predictions["model"].eq("physical_spatial_hgb")
    ].copy()
    style()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)

    ax = axes[0, 0]
    scatter = ax.scatter(frame["lon"], frame["lat"], c=frame["fold"], cmap="tab10", s=7, linewidths=0)
    ax.set_xlim(128, 146)
    ax.set_ylim(30, 46)
    ax.set_xlabel("Longitude (degrees E)")
    ax.set_ylabel("Latitude (degrees N)")
    ax.set_title("a  Five spatial validation blocks")
    fig.colorbar(scatter, ax=ax, ticks=range(5), label="Fold", fraction=0.046, pad=0.03)

    ax = axes[0, 1]
    fold = metrics[metrics["period_s"].eq(3.0) & metrics["scope"].eq("fold")].copy()
    models = ["physical_hgb", "physical_spatial_hgb"]
    labels = ["Site", "Site + location"]
    colors = ["#2A6F97", "#AA3377"]
    x = np.arange(5)
    for index, (model, label, color) in enumerate(zip(models, labels, colors, strict=True)):
        sub = fold[fold["model"].eq(model)].sort_values("fold")
        ax.bar(x + (index - 0.5) * 0.36, sub["rmse_reduction_vs_zero_pct"], 0.36, color=color, label=label)
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_xticks(x, [str(i) for i in x])
    ax.set_xlabel("Held-out spatial fold")
    ax.set_ylabel("RMSE reduction (%)")
    ax.set_title("b  SA(3.0 s) fold performance")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    x_obs = frame["station_effect_log10"].to_numpy(float)
    y_pred = frame["centered_oof_prediction_log10"].to_numpy(float)
    limit = float(np.quantile(np.abs(np.concatenate([x_obs, y_pred])), 0.995))
    ax.hexbin(x_obs, y_pred, gridsize=45, mincnt=1, cmap="Blues", linewidths=0)
    ax.plot([-limit, limit], [-limit, limit], color="#AA3377", lw=1.0)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("Observed station term (log10)")
    ax.set_ylabel("Out-of-fold prediction (log10)")
    ax.set_title(f"c  Held-out predictions (r = {np.corrcoef(x_obs, y_pred)[0, 1]:.3f})")

    ax = axes[1, 1]
    errors = [
        (frame.loc[frame["fold"].eq(fold_id), "centered_oof_prediction_log10"] - frame.loc[frame["fold"].eq(fold_id), "station_effect_log10"]).to_numpy()
        for fold_id in range(5)
    ]
    ax.boxplot(errors, tick_labels=[str(i) for i in range(5)], showfliers=False)
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_xlabel("Held-out spatial fold")
    ax.set_ylabel("Prediction error (log10)")
    ax.set_title("d  Error distributions")
    save(fig, "supplementary_figure_spatial_validation")


def event_repeatability_figure() -> None:
    repeat = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_holdout_station_repeatability.csv")
    decomposition = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_station_decomposition_metrics.csv")
    fold = repeat[repeat["scope"].eq("fold")]
    mean = repeat[repeat["scope"].eq("fold_mean")].sort_values("period_s")
    style()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.9), constrained_layout=True)

    ax = axes[0, 0]
    for fold_id, sub in fold.groupby("fold"):
        ax.plot(sub["period_s"], sub["train_test_station_correlation"], color="0.75", lw=0.8)
    ax.plot(mean["period_s"], mean["train_test_station_correlation"], color="#2A6F97", marker="o", lw=1.5)
    period_axis(ax)
    ax.set_ylim(0.85, 1.0)
    ax.set_ylabel("Train-test station correlation")
    ax.set_title("a  Station-term repeatability")

    ax = axes[0, 1]
    for fold_id, sub in fold.groupby("fold"):
        ax.plot(sub["period_s"], sub["rmse_reduction_vs_zero_pct"], color="0.75", lw=0.8)
    ax.plot(mean["period_s"], mean["rmse_reduction_vs_zero_pct"], color="#AA3377", marker="o", lw=1.5)
    period_axis(ax)
    ax.set_ylabel("Test RMSE reduction (%)")
    ax.set_title("b  Prediction across held-out events")

    ax = axes[1, 0]
    ax.plot(mean["period_s"], mean["n_paired_stations"], color="#2A6F97", marker="o", lw=1.4)
    period_axis(ax)
    ax.set_ylabel("Mean paired stations per fold")
    ax.set_title("c  Event-fold station support")

    ax = axes[1, 1]
    primary = decomposition[decomposition["model"].eq("mf2013_site")].drop_duplicates("period_s").sort_values("period_s")
    for column, label, color in [
        ("event_effect_weighted_std_log10", "Event", "#AA3377"),
        ("station_effect_weighted_std_log10", "Station", "#2A6F97"),
        ("remainder_rmse_log10", "Record remainder", "0.35"),
    ]:
        ax.plot(primary["period_s"], primary[column], marker="o", lw=1.3, color=color, label=label)
    period_axis(ax)
    ax.set_ylabel("Residual component (log10)")
    ax.set_title("d  Decomposition amplitudes")
    ax.legend(frameon=False)
    save(fig, "supplementary_figure_event_repeatability")


def response_spectrum_figure() -> None:
    values = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_adjusted_surface_spectrum_values.csv")
    summary = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_adjusted_surface_spectrum_summary.csv")
    predictions = pd.read_csv(SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv")
    preferred = predictions[predictions["model"].eq("physical_spatial_hgb")]
    style()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)

    ax = axes[0, 0]
    data = [preferred[preferred["period_s"].eq(period)]["oof_multiplier"].to_numpy() for period in PERIOD_TICKS]
    ax.boxplot(data, tick_labels=["0.1", "0.2", "0.3", "0.5", "1", "2", "3", "5"], showfliers=False)
    ax.axhline(1.0, color="0.3", lw=0.8)
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Cross-validated multiplier")
    ax.set_title("a  Station multiplier distributions")

    ax = axes[0, 1]
    rows = []
    for period, sub in preferred.groupby("period_s"):
        rows.append({"period_s": period, "q05": sub["oof_multiplier"].quantile(0.05), "q50": sub["oof_multiplier"].quantile(0.5), "q95": sub["oof_multiplier"].quantile(0.95)})
    quantiles = pd.DataFrame(rows).sort_values("period_s")
    ax.plot(quantiles["period_s"], quantiles["q50"], color="#2A6F97", marker="o", lw=1.4)
    ax.fill_between(quantiles["period_s"], quantiles["q05"], quantiles["q95"], color="#2A6F97", alpha=0.2)
    ax.axhline(1.0, color="0.3", lw=0.8)
    period_axis(ax)
    ax.set_ylabel("Multiplier")
    ax.set_title("b  Median and central 90% range")

    ax = axes[1, 0]
    sa3 = summary[summary["period_s"].eq(3.0)].copy()
    order = ["50y_39pct", "50y_10pct", "50y_5pct", "50y_2pct"]
    sa3["probability_level"] = pd.Categorical(sa3["probability_level"], order, ordered=True)
    sa3 = sa3.sort_values("probability_level")
    x = np.arange(len(sa3))
    for index, (column, label, color) in enumerate([
        ("official_vs400_sa_g_q50", "Official Vs400", "0.35"),
        ("ergodic_surface_sa_g_q50", "MF2013 surface", "#2A6F97"),
        ("adjusted_surface_sa_g_q50", "Station adjusted", "#AA3377"),
    ]):
        ax.plot(x, sa3[column], marker="o", lw=1.3, color=color, label=label)
    ax.set_xticks(x, ["39%", "10%", "5%", "2%"])
    ax.set_xlabel("50-year exceedance probability")
    ax.set_ylabel("Median SA(3.0 s) (g)")
    ax.set_title("c  Probability levels")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    matched = values[values["period_s"].eq(3.0) & values["probability_level"].eq("50y_10pct")]
    points = ax.scatter(
        matched["ergodic_surface_sa_g"],
        matched["adjusted_surface_sa_g"],
        c=matched["oof_multiplier"],
        cmap="RdBu_r",
        s=8,
        linewidths=0,
        alpha=0.65,
        rasterized=True,
    )
    low = float(min(matched["ergodic_surface_sa_g"].min(), matched["adjusted_surface_sa_g"].min()))
    high = float(max(matched["ergodic_surface_sa_g"].max(), matched["adjusted_surface_sa_g"].max()))
    ax.plot([low, high], [low, high], color="0.3", lw=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("MF2013 surface SA(3.0 s) (g)")
    ax.set_ylabel("Station-adjusted SA(3.0 s) (g)")
    ax.set_title("d  Matched-station ordinates")
    fig.colorbar(points, ax=ax, label="Multiplier", fraction=0.046, pad=0.03)
    save(fig, "supplementary_figure_response_spectra")


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    spatial_validation_figure()
    event_repeatability_figure()
    response_spectrum_figure()


if __name__ == "__main__":
    main()
