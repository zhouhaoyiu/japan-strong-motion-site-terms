#!/usr/bin/env python3
"""Estimate source-conditioned station fields and propagate them to J-SHIS maps."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core


ROOT = Path(__file__).resolve().parents[1]
PATH_TERMS = core.SUPPLEMENT_DIR / "jshis_path_stratification_station_terms.csv"
PRIMARY_TERMS = core.OUT_STATION_TERMS
ACTIVE_MAP = (
    ROOT
    / "work"
    / "external_data"
    / "jshis_respmap"
    / "P-Y2020-RESP-MAP-AVR-LND_MTTL-T50.zip"
)
SUBDUCTION_MAP = (
    ROOT
    / "work"
    / "external_data"
    / "jshis_respmap"
    / "P-Y2020-RESP-MAP-AVR-PPE_MTTL-T50.zip"
)

OUT_TERMS = core.SUPPLEMENT_DIR / "jshis_source_conditioned_station_terms.csv"
OUT_PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_source_conditioned_station_predictions.csv"
OUT_METRICS = core.SUPPLEMENT_DIR / "jshis_source_conditioned_station_metrics.csv"
OUT_BOOTSTRAP = core.SUPPLEMENT_DIR / "jshis_source_conditioned_station_bootstrap.csv"
OUT_HAZARD_VALUES = core.SUPPLEMENT_DIR / "jshis_source_conditioned_hazard_values.csv"
OUT_HAZARD_SUMMARY = core.SUPPLEMENT_DIR / "jshis_source_conditioned_hazard_summary.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_source_conditioned_station_fields.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "supplementary_figure_source_conditioned_fields.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "supplementary_figure_source_conditioned_fields.png"

SOURCE_TARGETS = {
    "Crustal": "source_conditioned_crustal",
    "Interplate": "source_conditioned_interplate",
    "Intraplate": "source_conditioned_intraplate",
}

HAZARD_SCENARIOS = {
    "active_shallow_crustal": ("Crustal", "active_shallow", ACTIVE_MAP),
    "subduction_interplate": ("Interplate", "subduction", SUBDUCTION_MAP),
    "subduction_intraplate": ("Intraplate", "subduction", SUBDUCTION_MAP),
}


def build_source_terms(
    path_terms: pd.DataFrame,
    primary_terms: pd.DataFrame,
    min_records: int,
) -> pd.DataFrame:
    metadata = primary_terms[primary_terms["model"].eq("mf2013_site")].drop_duplicates(
        ["siteid2", "period_s"]
    )
    metadata = metadata.drop(
        columns=["n_records", "station_effect_log10", "model", "response_component"]
    )
    source = path_terms[path_terms["stratification"].eq("source_type")].copy()
    source = source.rename(
        columns={
            "stratum_n_records": "n_records",
            "aligned_stratum_station_effect_log10": "station_effect_log10",
        }
    )
    source = source.merge(
        metadata,
        on=["siteid2", "period_s", "period_code"],
        how="left",
        validate="many_to_one",
    )
    source["response_component"] = "RotD100"
    source["source_type"] = source["stratum"]
    source["model"] = source["source_type"].map(SOURCE_TARGETS)
    source["meshcode250_str"] = source["meshcode250_str"].map(core.normalize_mesh)
    source["source_target_centering_shift_log10"] = np.nan
    for (source_type, period_s), index in source.groupby(
        ["source_type", "period_s"], sort=True
    ).groups.items():
        selected = source.loc[index, "n_records"].ge(min_records)
        supported_index = np.asarray(index)[selected.to_numpy()]
        if len(supported_index) == 0:
            raise ValueError(f"no supported stations for {source_type} at {period_s:g} s")
        weights = source.loc[supported_index, "n_records"].to_numpy(float)
        values = source.loc[supported_index, "station_effect_log10"].to_numpy(float)
        shift = float(np.average(values, weights=weights))
        source.loc[index, "station_effect_log10"] -= shift
        source.loc[index, "source_target_centering_shift_log10"] = shift
    return source


def bootstrap_metrics(
    predictions: pd.DataFrame,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    rows = []
    for group_index, ((source_type, period_s, model), block) in enumerate(
        predictions.groupby(["source_type", "period_s", "model"], sort=True)
    ):
        y = block["station_effect_log10"].to_numpy(float)
        predicted = block["centered_oof_prediction_log10"].to_numpy(float)
        weights = block["n_records"].to_numpy(float)
        rng = np.random.default_rng(seed + group_index)
        correlations = np.empty(replicates, dtype=float)
        gains = np.empty(replicates, dtype=float)
        for replicate in range(replicates):
            sampled = rng.integers(0, len(block), len(block))
            sample_y = y[sampled]
            sample_predicted = predicted[sampled]
            sample_weights = weights[sampled]
            correlations[replicate] = np.corrcoef(sample_y, sample_predicted)[0, 1]
            baseline = core.weighted_rmse(
                sample_y, np.zeros_like(sample_y), sample_weights
            )
            prediction_rmse = core.weighted_rmse(
                sample_y, sample_predicted, sample_weights
            )
            gains[replicate] = 100.0 * (baseline - prediction_rmse) / baseline
        rows.append(
            {
                "source_type": source_type,
                "period_s": period_s,
                "model": model,
                "n_stations": len(block),
                "bootstrap_replicates": replicates,
                "pearson_ci_low": float(np.nanquantile(correlations, 0.025)),
                "pearson_ci_high": float(np.nanquantile(correlations, 0.975)),
                "rmse_gain_pct_ci_low": float(np.nanquantile(gains, 0.025)),
                "rmse_gain_pct_ci_high": float(np.nanquantile(gains, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def build_hazard_outputs(
    predictions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    coefficients = core.load_coefficients(core.DEFAULT_COEFFICIENTS)
    value_blocks = []
    summary_blocks = []
    for scenario, (source_type, product, response_map) in HAZARD_SCENARIOS.items():
        source_predictions = predictions[predictions["source_type"].eq(source_type)].copy()
        values, summary = core.build_hazard_values(
            source_predictions, coefficients, response_map
        )
        values["hazard_scenario"] = scenario
        values["hazard_product"] = product
        values["source_type"] = source_type
        summary["hazard_scenario"] = scenario
        summary["hazard_product"] = product
        summary["source_type"] = source_type
        value_blocks.append(values)
        summary_blocks.append(summary)
        print(
            f"hazard propagation: {scenario} rows={len(values):,}",
            flush=True,
        )
    return pd.concat(value_blocks, ignore_index=True), pd.concat(
        summary_blocks, ignore_index=True
    )


def save_figure(
    predictions: pd.DataFrame,
    metrics: pd.DataFrame,
    hazard_summary: pd.DataFrame,
) -> None:
    selected_metrics = metrics[
        metrics["scope"].eq("overall")
        & metrics["model"].eq("physical_spatial_hgb")
    ].copy()
    colors = {
        "Crustal": "#D55E00",
        "Interplate": "#0072B2",
        "Intraplate": "#009E73",
    }
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.7,
            "savefig.dpi": 300,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6), constrained_layout=True)
    for source_type, block in selected_metrics.groupby("source_type", sort=False):
        block = block.sort_values("period_s")
        axes[0, 0].semilogx(
            block["period_s"],
            block["rmse_reduction_vs_zero_pct"],
            marker="o",
            color=colors[source_type],
            label=source_type,
        )
        axes[0, 1].semilogx(
            block["period_s"],
            block["observed_predicted_correlation"],
            marker="o",
            color=colors[source_type],
            label=source_type,
        )
    axes[0, 0].axhline(0.0, color="#555555", linestyle="--", linewidth=0.8)
    axes[0, 0].set(
        title="a  Source-conditioned spatial prediction",
        xlabel="Period (s)",
        ylabel="RMSE gain versus zero (%)",
    )
    axes[0, 1].axhline(0.0, color="#555555", linestyle="--", linewidth=0.8)
    axes[0, 1].set(
        title="b  Observed-predicted agreement",
        xlabel="Period (s)",
        ylabel="Pearson correlation",
    )
    axes[0, 1].legend(frameon=False)

    multiplier = predictions[
        predictions["period_s"].eq(3.0)
        & predictions["model"].eq("physical_spatial_hgb")
    ]
    groups = [
        multiplier.loc[multiplier["source_type"].eq(source), "oof_multiplier"]
        for source in SOURCE_TARGETS
    ]
    boxes = axes[1, 0].boxplot(
        groups,
        tick_labels=list(SOURCE_TARGETS),
        patch_artist=True,
        showfliers=False,
    )
    for box, source_type in zip(boxes["boxes"], SOURCE_TARGETS, strict=True):
        box.set(facecolor=colors[source_type], alpha=0.55)
    axes[1, 0].axhline(1.0, color="#555555", linestyle="--", linewidth=0.8)
    axes[1, 0].set(
        title="c  SA(3.0 s) out-of-fold multipliers",
        ylabel="Station multiplier",
    )

    selected_hazard = hazard_summary[hazard_summary["probability_level"].eq("50y_10pct")]
    scenario_labels = {
        "active_shallow_crustal": "Active shallow / crustal",
        "subduction_interplate": "Subduction / interplate",
        "subduction_intraplate": "Subduction / intraplate",
    }
    scenario_colors = {
        "active_shallow_crustal": colors["Crustal"],
        "subduction_interplate": colors["Interplate"],
        "subduction_intraplate": colors["Intraplate"],
    }
    for scenario, block in selected_hazard.groupby("hazard_scenario", sort=False):
        block = block.sort_values("period_s")
        axes[1, 1].loglog(
            block["period_s"],
            block["adjusted_surface_sa_g_q50"],
            marker="o",
            color=scenario_colors[scenario],
            label=scenario_labels[scenario],
        )
    axes[1, 1].set(
        title="d  Source-conditioned surface spectra",
        xlabel="Period (s)",
        ylabel="Median adjusted surface SA (g)",
    )
    axes[1, 1].legend(frameon=False)
    for axis in axes.flat:
        axis.grid(alpha=0.2, linewidth=0.6)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    terms: pd.DataFrame,
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    hazard_summary: pd.DataFrame,
    min_records: int,
    spatial_folds: int,
) -> None:
    selected_metrics = metrics[
        metrics["scope"].eq("overall")
        & metrics["model"].eq("physical_spatial_hgb")
        & metrics["period_s"].isin([1.0, 2.0, 3.0, 5.0])
    ]
    selected_bootstrap = bootstrap[
        bootstrap["model"].eq("physical_spatial_hgb")
        & bootstrap["period_s"].isin([1.0, 2.0, 3.0, 5.0])
    ]
    lines = [
        "# Source-conditioned station-field and response-spectrum sensitivity",
        "",
        "## Locked design",
        "",
        "- Source classes are the public flatfile labels: crustal, interplate and intraplate. Every class, period and spatial fold is retained.",
        f"- A station requires at least {min_records} records within a source class. Targets are record-weighted to zero within source class and period before fitting.",
        f"- The fixed primary site-plus-location model and {spatial_folds} spatial blocks are used without tuning. Confidence intervals use 2,000 station-cluster bootstrap replicates.",
        "- The crustal field is applied to the official active-shallow product. Interplate and intraplate fields are applied separately to the official subduction product as two sensitivities; they are not averaged because official source weights are unavailable.",
        "",
        "## Spatial prediction",
        "",
    ]
    merged = selected_metrics.merge(
        selected_bootstrap,
        on=["source_type", "period_s", "model", "n_stations"],
        how="left",
        validate="one_to_one",
    )
    for source_type in SOURCE_TARGETS:
        lines.append(f"### {source_type}")
        lines.append("")
        for _, row in merged[merged["source_type"].eq(source_type)].sort_values(
            "period_s"
        ).iterrows():
            lines.append(
                f"- {row['period_s']:g} s: {int(row['n_stations']):,} stations; "
                f"correlation {row['observed_predicted_correlation']:.3f} "
                f"({row['pearson_ci_low']:.3f}--{row['pearson_ci_high']:.3f}); "
                f"RMSE gain {row['rmse_reduction_vs_zero_pct']:.1f}% "
                f"({row['rmse_gain_pct_ci_low']:.1f}%--{row['rmse_gain_pct_ci_high']:.1f}%)."
            )
        lines.append("")
    lines.extend(["## Response-spectrum propagation", ""])
    selected_hazard = hazard_summary[
        hazard_summary["period_s"].eq(3.0)
        & hazard_summary["probability_level"].eq("50y_10pct")
    ]
    for _, row in selected_hazard.sort_values("hazard_scenario").iterrows():
        lines.append(
            f"- {row['hazard_scenario']}: {int(row['n_station_values']):,} matched stations; "
            f"median official {row['official_vs400_sa_g_q50']:.3f} g; "
            f"median surface reference {row['ergodic_surface_sa_g_q50']:.3f} g; "
            f"median adjusted {row['adjusted_surface_sa_g_q50']:.3f} g."
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "These fields condition the empirical station residual on the observed source class. They retain unresolved path effects within each class and are sensitivities on published fixed-probability products. They do not replace source-specific terms inside a production PSHA integral.",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_outputs(
    terms: pd.DataFrame,
    predictions: pd.DataFrame,
    metrics: pd.DataFrame,
    bootstrap: pd.DataFrame,
    hazard_summary: pd.DataFrame,
) -> None:
    periods = {item.period_s for item in core.PERIODS}
    assert set(terms["source_type"]) == set(SOURCE_TARGETS)
    assert set(predictions["source_type"]) == set(SOURCE_TARGETS)
    assert set(metrics["source_type"]) == set(SOURCE_TARGETS)
    assert set(bootstrap["source_type"]) == set(SOURCE_TARGETS)
    assert set(predictions["period_s"]) == periods
    assert set(hazard_summary["hazard_scenario"]) == set(HAZARD_SCENARIOS)
    assert predictions["centered_oof_prediction_log10"].notna().all()
    assert bootstrap[
        [
            "pearson_ci_low",
            "pearson_ci_high",
            "rmse_gain_pct_ci_low",
            "rmse_gain_pct_ci_high",
        ]
    ].notna().all().all()
    for path in [
        OUT_TERMS,
        OUT_PREDICTIONS,
        OUT_METRICS,
        OUT_BOOTSTRAP,
        OUT_HAZARD_VALUES,
        OUT_HAZARD_SUMMARY,
        OUT_AUDIT,
        OUT_FIGURE_PDF,
        OUT_FIGURE_PNG,
    ]:
        assert path.exists() and path.stat().st_size > 0, path


def run(args: argparse.Namespace) -> None:
    core.SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    core.FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    path_terms = pd.read_csv(args.path_terms)
    primary_terms = pd.read_csv(args.primary_terms)
    terms = build_source_terms(path_terms, primary_terms, args.min_records)
    prediction_blocks = []
    metric_blocks = []
    for source_type, target in SOURCE_TARGETS.items():
        predictions, metrics = core.cross_validate_station_models(
            terms,
            min_records=args.min_records,
            n_splits=args.spatial_folds,
            seed=args.seed,
            target_model=target,
        )
        predictions["source_type"] = source_type
        metrics["source_type"] = source_type
        prediction_blocks.append(predictions)
        metric_blocks.append(metrics)
    predictions = pd.concat(prediction_blocks, ignore_index=True)
    metrics = pd.concat(metric_blocks, ignore_index=True)
    bootstrap = bootstrap_metrics(predictions, args.bootstrap_replicates, args.seed)
    terms.to_csv(OUT_TERMS, index=False)
    predictions.to_csv(OUT_PREDICTIONS, index=False)
    metrics.to_csv(OUT_METRICS, index=False)
    bootstrap.to_csv(OUT_BOOTSTRAP, index=False)

    if not predictions["meshcode250_str"].astype(str).str.fullmatch(r"\d{10}").all():
        raise ValueError("source-conditioned predictions contain invalid 250 m mesh codes")
    hazard_values, hazard_summary = build_hazard_outputs(predictions)
    hazard_values.to_csv(OUT_HAZARD_VALUES, index=False)
    hazard_summary.to_csv(OUT_HAZARD_SUMMARY, index=False)
    save_figure(predictions, metrics, hazard_summary)
    write_audit(
        terms,
        metrics,
        bootstrap,
        hazard_summary,
        min_records=args.min_records,
        spatial_folds=args.spatial_folds,
    )
    validate_outputs(terms, predictions, metrics, bootstrap, hazard_summary)
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-terms", type=Path, default=PATH_TERMS)
    parser.add_argument("--primary-terms", type=Path, default=PRIMARY_TERMS)
    parser.add_argument("--min-records", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=20260711)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
