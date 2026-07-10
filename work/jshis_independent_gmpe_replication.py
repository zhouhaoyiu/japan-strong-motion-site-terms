#!/usr/bin/env python3
"""Independent Zhao et al. (2006) replication of the station-term result.

The script applies the OpenQuake hazardlib implementation of Zhao et al. (2006)
to the same J-SHIS/NIED records and eight direct response-spectrum periods used
by the primary MF2013 analysis. It then repeats the event-station decomposition,
event-group holdout, and spatial-block station prediction without changing the
validation design.

The flatfile's documented shortest fault distance is mapped to Zhao's rupture
distance input, and the public AVS30 field supplies Vs30. This is an independent
ground-motion-model sensitivity test, not an independent observational dataset.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import math
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openquake.hazardlib.gsim.zhao_2006 import (
    ZhaoEtAl2006Asc,
    ZhaoEtAl2006SInter,
    ZhaoEtAl2006SSlab,
)
from openquake.hazardlib.imt import SA

import jshis_event_adjusted_station_model as core


ARTICLE_DIR = core.ARTICLE_DIR
SUPPLEMENT_DIR = core.SUPPLEMENT_DIR
FIGURE_DIR = core.FIGURE_DIR

OUT_STATION_TERMS = SUPPLEMENT_DIR / "jshis_zhao2006_station_terms.csv"
OUT_DECOMPOSITION = SUPPLEMENT_DIR / "jshis_zhao2006_decomposition_metrics.csv"
OUT_EVENT_HOLDOUT = SUPPLEMENT_DIR / "jshis_zhao2006_event_holdout_repeatability.csv"
OUT_MODEL_PREDICTIONS = SUPPLEMENT_DIR / "jshis_zhao2006_station_model_predictions.csv"
OUT_MODEL_METRICS = SUPPLEMENT_DIR / "jshis_zhao2006_station_model_metrics.csv"
OUT_REPLICATION = SUPPLEMENT_DIR / "jshis_independent_gmpe_station_term_replication.csv"
OUT_AUDIT = SUPPLEMENT_DIR / "jshis_independent_gmpe_replication.md"
OUT_FIGURE_PDF = FIGURE_DIR / "figure_ground_motion_model_sensitivity.pdf"
OUT_FIGURE_PNG = FIGURE_DIR / "figure_ground_motion_model_sensitivity.png"

MF2013_TERMS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_terms.csv"
MF2013_EVENT_HOLDOUT = SUPPLEMENT_DIR / "jshis_event_holdout_station_repeatability.csv"
MF2013_MODEL_METRICS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_metrics.csv"

SOURCE_MODELS = {
    1: ZhaoEtAl2006Asc(),
    2: ZhaoEtAl2006SInter(),
    3: ZhaoEtAl2006SSlab(),
}


def build_context(frame: pd.DataFrame) -> np.recarray:
    context = np.recarray(
        len(frame),
        dtype=[
            ("mag", float),
            ("hypo_depth", float),
            ("rrup", float),
            ("vs30", float),
            ("rake", float),
        ],
    )
    context.mag = pd.to_numeric(frame["mw"], errors="coerce").to_numpy(float)
    context.hypo_depth = pd.to_numeric(frame["jem_depth"], errors="coerce").to_numpy(float)
    context.rrup = pd.to_numeric(frame["fault_dist"], errors="coerce").clip(lower=0.1).to_numpy(float)
    context.vs30 = pd.to_numeric(frame["avs30"], errors="coerce").to_numpy(float)
    context.rake = pd.to_numeric(frame["rake1"], errors="coerce").fillna(0.0).to_numpy(float)
    return context


def predict_chunk(frame: pd.DataFrame) -> np.ndarray:
    imts = [SA(spec.period_s) for spec in core.PERIODS]
    predictions = np.full((len(core.PERIODS), len(frame)), np.nan, dtype=float)
    for source_class, model in SOURCE_MODELS.items():
        mask = frame["eq_location_type_id"].eq(source_class)
        if source_class == 1:
            mask &= pd.to_numeric(frame["rake1"], errors="coerce").notna()
        indices = np.flatnonzero(mask.to_numpy())
        if not len(indices):
            continue
        context = build_context(frame.iloc[indices])
        mean = np.zeros((len(imts), len(indices)), dtype=float)
        sigma = np.zeros_like(mean)
        tau = np.zeros_like(mean)
        phi = np.zeros_like(mean)
        model.compute(context, imts, mean, sigma, tau, phi)
        predictions[:, indices] = mean / math.log(10.0) + math.log10(core.G_IN_CM_S2)
    return predictions


def load_zhao_records(
    flatfile: Path,
    site: pd.DataFrame,
    source: pd.DataFrame,
    chunksize: int,
) -> dict[float, pd.DataFrame]:
    usecols = ["smrec_id", "site_id", "eq_source_id", "fault_dist", *[spec.rotd50_col for spec in core.PERIODS]]
    buffers: dict[float, list[pd.DataFrame]] = {spec.period_s: [] for spec in core.PERIODS}
    joined_rows = 0
    with zipfile.ZipFile(flatfile) as archive:
        reader = pd.read_csv(
            archive.open("smrec_schema.tsv"),
            sep="\t",
            usecols=usecols,
            chunksize=chunksize,
            low_memory=False,
        )
        for chunk_number, chunk in enumerate(reader, start=1):
            chunk["site_id"] = pd.to_numeric(chunk["site_id"], errors="coerce").astype("Int64")
            chunk["siteid2"] = (chunk["site_id"] // 10).astype("Int64")
            chunk["eq_source_id"] = pd.to_numeric(chunk["eq_source_id"], errors="coerce").astype("Int64")
            joined = chunk.merge(site, on="siteid2", how="left", validate="many_to_one")
            joined = joined.merge(source, on="eq_source_id", how="left", validate="many_to_one")
            prediction = predict_chunk(joined)
            joined_rows += len(joined)
            common = (
                pd.to_numeric(joined["fault_dist"], errors="coerce").gt(0)
                & pd.to_numeric(joined["avs30"], errors="coerce").gt(0)
                & pd.to_numeric(joined["mw"], errors="coerce").notna()
                & pd.to_numeric(joined["jem_depth"], errors="coerce").notna()
                & joined["eq_location_type_id"].isin([1, 2, 3])
                & joined["installation_situation_id"].eq(core.PRIMARY_INSTALLATION_SITUATION_ID)
            )
            common &= ~(
                joined["eq_location_type_id"].eq(1)
                & pd.to_numeric(joined["rake1"], errors="coerce").isna()
            )
            for period_index, spec in enumerate(core.PERIODS):
                observed_raw = pd.to_numeric(joined[spec.rotd50_col], errors="coerce")
                observed = np.log10(observed_raw.where(observed_raw > 0))
                predicted = pd.Series(prediction[period_index], index=joined.index)
                usable = common & observed.notna() & predicted.notna() & np.isfinite(observed) & np.isfinite(predicted)
                block = joined.loc[usable, ["smrec_id", "siteid2", "eq_source_id"]].copy()
                block["residual_site"] = observed.loc[usable].to_numpy() - predicted.loc[usable].to_numpy()
                buffers[spec.period_s].append(block)
            print(f"Zhao records: chunk={chunk_number} joined={joined_rows:,}", flush=True)
    return {period: pd.concat(parts, ignore_index=True) for period, parts in buffers.items()}


def decompose_zhao(
    records: dict[float, pd.DataFrame],
    site: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    term_rows = []
    metric_rows = []
    for spec in core.PERIODS:
        terms, _, metrics = core.two_way_decomposition(records[spec.period_s], "residual_site")
        terms["period_s"] = spec.period_s
        terms["period_code"] = spec.code
        terms["model"] = "zhao2006_openquake"
        terms["response_component"] = "RotD50"
        term_rows.append(terms)
        metric_rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "model": "zhao2006_openquake",
                "response_component": "RotD50",
                **metrics,
            }
        )
        print(f"Zhao decomposition: period={spec.period_s:g}s", flush=True)
    station_terms = pd.concat(term_rows, ignore_index=True).merge(
        site.drop_duplicates("siteid2"),
        on="siteid2",
        how="left",
        validate="many_to_one",
    )
    return station_terms, pd.DataFrame(metric_rows)


def weighted_correlation(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    x_mean = np.average(x, weights=weights)
    y_mean = np.average(y, weights=weights)
    covariance = np.average((x - x_mean) * (y - y_mean), weights=weights)
    x_var = np.average((x - x_mean) ** 2, weights=weights)
    y_var = np.average((y - y_mean) ** 2, weights=weights)
    return float(covariance / np.sqrt(x_var * y_var))


def compare_station_terms(zhao_terms: pd.DataFrame, min_records: int) -> pd.DataFrame:
    mf2013 = pd.read_csv(MF2013_TERMS)
    mf2013 = mf2013[
        mf2013["model"].eq("mf2013_site") & mf2013["n_records"].ge(min_records)
    ][["siteid2", "period_s", "n_records", "station_effect_log10"]].rename(
        columns={
            "n_records": "mf2013_n_records",
            "station_effect_log10": "mf2013_station_effect_log10",
        }
    )
    zhao = zhao_terms[zhao_terms["n_records"].ge(min_records)][
        ["siteid2", "period_s", "n_records", "station_effect_log10"]
    ].rename(
        columns={
            "n_records": "zhao_n_records",
            "station_effect_log10": "zhao_station_effect_log10",
        }
    )
    paired = mf2013.merge(zhao, on=["siteid2", "period_s"], how="inner", validate="one_to_one")
    rows = []
    for spec in core.PERIODS:
        frame = paired[paired["period_s"].eq(spec.period_s)].copy()
        x = frame["mf2013_station_effect_log10"].to_numpy(float)
        y = frame["zhao_station_effect_log10"].to_numpy(float)
        weights = np.minimum(frame["mf2013_n_records"], frame["zhao_n_records"]).to_numpy(float)
        threshold = np.quantile(np.abs(x), 0.5)
        strong = (np.abs(x) >= threshold) & (np.abs(y) >= threshold)
        rows.append(
            {
                "period_s": spec.period_s,
                "period_code": spec.code,
                "n_paired_stations": len(frame),
                "pearson_correlation": float(np.corrcoef(x, y)[0, 1]),
                "weighted_pearson_correlation": weighted_correlation(x, y, weights),
                "spearman_correlation": float(pd.Series(x).corr(pd.Series(y), method="spearman")),
                "median_absolute_difference_log10": float(np.median(np.abs(x - y))),
                "sign_agreement_pct": float(100.0 * np.mean(np.sign(x) == np.sign(y))),
                "strong_term_sign_agreement_pct": float(100.0 * np.mean(np.sign(x[strong]) == np.sign(y[strong]))),
            }
        )
    return pd.DataFrame(rows)


def save_figure(
    zhao_terms: pd.DataFrame,
    replication: pd.DataFrame,
    zhao_event: pd.DataFrame,
    zhao_metrics: pd.DataFrame,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    mf_terms = pd.read_csv(MF2013_TERMS)
    mf3 = mf_terms[
        mf_terms["model"].eq("mf2013_site")
        & mf_terms["period_s"].eq(3.0)
        & mf_terms["n_records"].ge(20)
    ][["siteid2", "station_effect_log10"]].rename(columns={"station_effect_log10": "mf2013"})
    zhao3 = zhao_terms[
        zhao_terms["period_s"].eq(3.0) & zhao_terms["n_records"].ge(20)
    ][["siteid2", "station_effect_log10"]].rename(columns={"station_effect_log10": "zhao"})
    paired = mf3.merge(zhao3, on="siteid2", validate="one_to_one")

    mf_event = pd.read_csv(MF2013_EVENT_HOLDOUT)
    mf_event = mf_event[mf_event["scope"].eq("fold_mean")]
    z_event = zhao_event[zhao_event["scope"].eq("fold_mean")]
    mf_model = pd.read_csv(MF2013_MODEL_METRICS)
    mf_model = mf_model[
        mf_model["scope"].eq("overall") & mf_model["model"].eq("physical_spatial_hgb")
    ]
    z_model = zhao_metrics[
        zhao_metrics["scope"].eq("overall") & zhao_metrics["model"].eq("physical_spatial_hgb")
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.8), constrained_layout=True)
    ax = axes[0, 0]
    ax.scatter(paired["mf2013"], paired["zhao"], s=8, alpha=0.25, color="#277DA1", edgecolors="none")
    limits = np.quantile(np.abs(np.r_[paired["mf2013"], paired["zhao"]]), 0.995)
    ax.plot([-limits, limits], [-limits, limits], color="#E63946", lw=1.2)
    ax.set(xlim=(-limits, limits), ylim=(-limits, limits), xlabel="MF2013 station term (log10)", ylabel="Zhao 2006 station term (log10)")
    r3 = replication.loc[replication["period_s"].eq(3.0), "pearson_correlation"].iloc[0]
    ax.set_title(f"a  SA(3.0 s) station terms (r = {r3:.3f})", loc="left")

    ax = axes[0, 1]
    ax.plot(replication["period_s"], replication["pearson_correlation"], "o-", color="#277DA1", label="Pearson")
    ax.plot(replication["period_s"], replication["spearman_correlation"], "s-", color="#E76F51", label="Spearman")
    ax.set_xscale("log")
    ax.set_xticks([spec.period_s for spec in core.PERIODS], [f"{spec.period_s:g}" for spec in core.PERIODS])
    ax.set(xlabel="Period (s)", ylabel="MF2013-Zhao station-term correlation", ylim=(0, 1.02))
    ax.legend(frameon=False)
    ax.set_title("b  Independent-GMPE agreement", loc="left")

    ax = axes[1, 0]
    ax.plot(mf_event["period_s"], mf_event["train_test_station_correlation"], "o-", color="#277DA1", label="MF2013")
    ax.plot(z_event["period_s"], z_event["train_test_station_correlation"], "s-", color="#E76F51", label="Zhao 2006")
    ax.set_xscale("log")
    ax.set_xticks([spec.period_s for spec in core.PERIODS], [f"{spec.period_s:g}" for spec in core.PERIODS])
    ax.set(xlabel="Period (s)", ylabel="Event-holdout station correlation", ylim=(0.85, 1.0))
    ax.legend(frameon=False)
    ax.set_title("c  Repeatability across held-out events", loc="left")

    ax = axes[1, 1]
    ax.plot(mf_model["period_s"], mf_model["rmse_reduction_vs_zero_pct"], "o-", color="#277DA1", label="MF2013")
    ax.plot(z_model["period_s"], z_model["rmse_reduction_vs_zero_pct"], "s-", color="#E76F51", label="Zhao 2006")
    ax.axhline(0, color="#555555", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks([spec.period_s for spec in core.PERIODS], [f"{spec.period_s:g}" for spec in core.PERIODS])
    ax.set(xlabel="Period (s)", ylabel="Spatial-block RMSE reduction (%)")
    ax.legend(frameon=False)
    ax.set_title("d  Public-variable spatial prediction", loc="left")

    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    records: dict[float, pd.DataFrame],
    decomposition: pd.DataFrame,
    replication: pd.DataFrame,
    event_holdout: pd.DataFrame,
    model_metrics: pd.DataFrame,
) -> None:
    sa3_replication = replication[replication["period_s"].eq(3.0)].iloc[0]
    sa3_event = event_holdout[
        event_holdout["scope"].eq("fold_mean") & event_holdout["period_s"].eq(3.0)
    ].iloc[0]
    sa3_model = model_metrics[
        model_metrics["scope"].eq("overall")
        & model_metrics["model"].eq("physical_spatial_hgb")
        & model_metrics["period_s"].eq(3.0)
    ].iloc[0]
    lines = [
        "# Independent GMPE replication",
        "",
        f"OpenQuake engine version: {importlib.metadata.version('openquake.engine')}.",
        "",
        "## Input mapping",
        "",
        "- Zhao et al. (2006) active-crust, subduction-interface and subduction-slab variants are selected from the J-SHIS source class.",
        "- The documented J-SHIS shortest fault distance is supplied as rupture distance.",
        "- J-SHIS AVS30 is supplied as Vs30 for every station.",
        "- Observed RotD50 spectral acceleration and Zhao predictions are compared in log10(cm/s2).",
        "",
        "## Coverage",
        "",
        *[f"- {spec.period_s:g} s: {len(records[spec.period_s]):,} records." for spec in core.PERIODS],
        "",
        "## SA(3.0 s) result",
        "",
        f"- MF2013-Zhao station-term Pearson correlation: {sa3_replication['pearson_correlation']:.3f}.",
        f"- MF2013-Zhao station-term Spearman correlation: {sa3_replication['spearman_correlation']:.3f}.",
        f"- Event-holdout station correlation: {sa3_event['train_test_station_correlation']:.3f}.",
        f"- Spatial-block RMSE reduction: {sa3_model['rmse_reduction_vs_zero_pct']:.1f}%.",
        "",
        "## Boundary",
        "",
        "This analysis changes the regional ground-motion model while retaining the same observations. It tests model dependence of the station field; it is not an independent-network validation or an official Zhao input-preparation reproduction.",
        "",
        "## Convergence",
        "",
        f"- Maximum final fixed-effect change: {decomposition['max_parameter_change'].max():.3e} log10 units.",
        f"- Maximum iteration count: {int(decomposition['iterations'].max())}.",
    ]
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    site, source = core.load_metadata(args.flatfile)
    records = load_zhao_records(args.flatfile, site, source, args.chunksize)
    station_terms, decomposition = decompose_zhao(records, site)
    replication = compare_station_terms(station_terms, args.min_station_records)
    event_holdout = core.event_holdout_repeatability(records, n_splits=args.spatial_folds, seed=args.seed)
    event_holdout["residual_model"] = "zhao2006_openquake"
    predictions, model_metrics = core.cross_validate_station_models(
        station_terms,
        min_records=args.min_station_records,
        n_splits=args.spatial_folds,
        seed=args.seed,
        target_model="zhao2006_openquake",
    )
    station_terms.to_csv(OUT_STATION_TERMS, index=False)
    decomposition.to_csv(OUT_DECOMPOSITION, index=False)
    event_holdout.to_csv(OUT_EVENT_HOLDOUT, index=False)
    predictions.to_csv(OUT_MODEL_PREDICTIONS, index=False)
    model_metrics.to_csv(OUT_MODEL_METRICS, index=False)
    replication.to_csv(OUT_REPLICATION, index=False)
    save_figure(station_terms, replication, event_holdout, model_metrics)
    write_audit(records, decomposition, replication, event_holdout, model_metrics)
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flatfile", type=Path, default=core.DEFAULT_FLATFILE)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--min-station-records", type=int, default=20)
    parser.add_argument("--spatial-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
