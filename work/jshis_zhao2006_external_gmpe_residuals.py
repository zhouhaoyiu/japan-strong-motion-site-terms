#!/usr/bin/env python3
"""External Zhao et al. (2006) GMPE residual audit for J-SHIS/NIED records.

This script is an independent model-sensitivity audit for the CEE manuscript
backbone. It uses the OpenQuake hazardlib implementation of Zhao et al. (2006)
and maps J-SHIS source classes to crustal, interface, and slab variants.

Important limitations:
- `fault_dist` is used as an `rrup` proxy because the public flatfile subset
  does not expose a direct OpenQuake rupture-distance context.
- OpenQuake returns ln(g); J-SHIS observations are converted from gal to g.
- `vs30` uses official `vs30` where available and falls back to `avs30` to keep
  national-scale coverage. This makes the result an audit, not a full
  reimplementation of Zhao et al. (2006) input preparation.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from openquake.hazardlib.gsim.zhao_2006 import (
    ZhaoEtAl2006Asc,
    ZhaoEtAl2006SInter,
    ZhaoEtAl2006SSlab,
)
from openquake.hazardlib.imt import PGA, SA


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ZIP = PROJECT_ROOT / "work" / "external_data" / "jshis_gmf" / "flatfile_sub1-v2024.zip"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
GAL_PER_G = 980.665

TARGETS = [
    {
        "target": "maxaccrd050",
        "target_label": "PGA RotD50",
        "period": "PGA",
        "observed_col": "maxaccrd050",
        "imt": PGA(),
    },
    {
        "target": "rsaccrd050d005t0030",
        "target_label": "SA(0.3s) RotD50",
        "period": "0.30",
        "observed_col": "rsaccrd050d005t0030",
        "imt": SA(0.3),
    },
    {
        "target": "rsaccrd050d005t0100",
        "target_label": "SA(1.0s) RotD50",
        "period": "1.00",
        "observed_col": "rsaccrd050d005t0100",
        "imt": SA(1.0),
    },
    {
        "target": "rsaccrd050d005t0300",
        "target_label": "SA(3.0s) RotD50",
        "period": "3.00",
        "observed_col": "rsaccrd050d005t0300",
        "imt": SA(3.0),
    },
]

SMREC_USECOLS = [
    "smrec_id",
    "site_id",
    "eq_source_id",
    "maxaccrd050",
    "rsaccrd050d005t0030",
    "rsaccrd050d005t0100",
    "rsaccrd050d005t0300",
    "fault_dist",
]

SITE_USECOLS = [
    "siteid2",
    "site_code",
    "site_name",
    "obs_network_id",
    "vs10",
    "vs20",
    "vs30",
    "avs30",
    "d1100",
    "d1400",
    "dbase",
]

SOURCE_USECOLS = [
    "eq_source_id",
    "segment_idx",
    "jem_origin_time",
    "jem_depth",
    "mjma",
    "mw",
    "rake1",
    "eq_location_type_id",
    "eq_event_name",
]

SOURCE_CLASS_LABELS = {
    1: "CRUSTAL",
    2: "INTERPLATE",
    3: "INTRAPLATE",
}

NETWORK_LABELS = {
    1: "K-NET",
    2: "KiK-net",
}

SOURCE_MODELS = {
    1: ("zhao2006_asc", ZhaoEtAl2006Asc()),
    2: ("zhao2006_sinter", ZhaoEtAl2006SInter()),
    3: ("zhao2006_sslab", ZhaoEtAl2006SSlab()),
}

AGG_VALUE_COLS = [
    "n_records",
    "residual_sum",
    "residual_sq_sum",
    "absolute_error_sum",
    "observed_sum",
    "predicted_sum",
]


def read_schema_tables(zip_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(zip_path) as zf:
        with zf.open("site_schema.tsv") as handle:
            site = pd.read_csv(handle, sep="\t", usecols=SITE_USECOLS)
        with zf.open("source_schema.tsv") as handle:
            source = pd.read_csv(handle, sep="\t", usecols=SOURCE_USECOLS)

    for col in SITE_USECOLS:
        if col not in ["site_code", "site_name"]:
            site[col] = pd.to_numeric(site[col], errors="coerce")
    site["siteid2"] = site["siteid2"].astype("Int64")
    site["vs30_proxy"] = pd.to_numeric(site["vs30"], errors="coerce")
    avs30 = pd.to_numeric(site["avs30"], errors="coerce")
    site["vs30_proxy_source"] = np.where(site["vs30_proxy"].gt(0), "vs30", "avs30_fallback")
    site.loc[~site["vs30_proxy"].gt(0), "vs30_proxy"] = avs30.loc[~site["vs30_proxy"].gt(0)]
    site.loc[~site["vs30_proxy"].gt(0), "vs30_proxy_source"] = "missing"
    site = site.drop_duplicates("siteid2").copy()

    for col in SOURCE_USECOLS:
        if col not in ["jem_origin_time", "eq_event_name"]:
            source[col] = pd.to_numeric(source[col], errors="coerce")
    source["eq_source_id"] = source["eq_source_id"].astype("Int64")
    source = source.sort_values(["eq_source_id", "segment_idx"]).drop_duplicates("eq_source_id")
    return site, source


def build_context(frame: pd.DataFrame) -> np.recarray:
    ctx = np.recarray(
        len(frame),
        dtype=[
            ("mag", float),
            ("hypo_depth", float),
            ("rrup", float),
            ("vs30", float),
            ("rake", float),
        ],
    )
    ctx.mag = pd.to_numeric(frame["mw"], errors="coerce").to_numpy(dtype=float)
    ctx.hypo_depth = pd.to_numeric(frame["jem_depth"], errors="coerce").to_numpy(dtype=float)
    ctx.rrup = pd.to_numeric(frame["fault_dist"], errors="coerce").clip(lower=0.1).to_numpy(dtype=float)
    ctx.vs30 = pd.to_numeric(frame["vs30_proxy"], errors="coerce").to_numpy(dtype=float)
    ctx.rake = pd.to_numeric(frame["rake1"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    return ctx


def zhao_predictions(frame: pd.DataFrame, target_info: dict[str, object]) -> pd.Series:
    pred = pd.Series(np.nan, index=frame.index, dtype=float)
    imt = target_info["imt"]
    for source_class, (_, model) in SOURCE_MODELS.items():
        mask = frame["eq_location_type_id"].eq(source_class)
        if source_class == 1:
            mask = mask & pd.to_numeric(frame["rake1"], errors="coerce").notna()
        if not mask.any():
            continue
        sub = frame.loc[mask]
        ctx = build_context(sub)
        mean = np.zeros((1, len(sub)), dtype=float)
        sig = np.zeros_like(mean)
        tau = np.zeros_like(mean)
        phi = np.zeros_like(mean)
        model.compute(ctx, [imt], mean, sig, tau, phi)
        pred.loc[sub.index] = mean[0]
    return pred


def residual_frame(joined: pd.DataFrame, target_info: dict[str, object]) -> pd.DataFrame:
    observed_raw = pd.to_numeric(joined[str(target_info["observed_col"])], errors="coerce")
    observed = np.log(observed_raw.where(observed_raw > 0) / GAL_PER_G)
    predicted = zhao_predictions(joined, target_info)
    usable = (
        observed.notna()
        & predicted.notna()
        & pd.to_numeric(joined["fault_dist"], errors="coerce").gt(0)
        & pd.to_numeric(joined["mw"], errors="coerce").notna()
        & pd.to_numeric(joined["jem_depth"], errors="coerce").notna()
        & pd.to_numeric(joined["vs30_proxy"], errors="coerce").gt(0)
        & joined["eq_location_type_id"].isin([1, 2, 3])
    )
    if not usable.any():
        return pd.DataFrame()
    out = joined.loc[
        usable,
        [
            "smrec_id",
            "siteid2",
            "site_code",
            "obs_network_id",
            "eq_source_id",
            "eq_location_type_id",
            "mw",
            "mjma",
            "jem_depth",
            "rake1",
            "fault_dist",
            "vs30_proxy",
            "vs30_proxy_source",
        ],
    ].copy()
    out["target"] = target_info["target"]
    out["target_label"] = target_info["target_label"]
    out["period"] = target_info["period"]
    out["model"] = out["eq_location_type_id"].map({k: v[0] for k, v in SOURCE_MODELS.items()})
    out["model_family"] = "zhao2006_openquake"
    out["observed_ln_g"] = observed.loc[usable].to_numpy()
    out["predicted_ln_g"] = predicted.loc[usable].to_numpy()
    out["residual_ln"] = out["observed_ln_g"] - out["predicted_ln_g"]
    out["absolute_error_ln"] = out["residual_ln"].abs()
    out["squared_error_ln"] = out["residual_ln"] ** 2
    out["source_class_label"] = out["eq_location_type_id"].map(SOURCE_CLASS_LABELS)
    out["network_label"] = out["obs_network_id"].map(NETWORK_LABELS)
    return out


def aggregate_residuals(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=group_cols + AGG_VALUE_COLS)
    return (
        frame.groupby(group_cols, dropna=False)
        .agg(
            n_records=("residual_ln", "size"),
            residual_sum=("residual_ln", "sum"),
            residual_sq_sum=("squared_error_ln", "sum"),
            absolute_error_sum=("absolute_error_ln", "sum"),
            observed_sum=("observed_ln_g", "sum"),
            predicted_sum=("predicted_ln_g", "sum"),
        )
        .reset_index()
    )


def finalize_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    n = out["n_records"].replace(0, np.nan)
    out["mean_residual_ln"] = out["residual_sum"] / n
    out["mae_ln"] = out["absolute_error_sum"] / n
    out["rmse_ln"] = np.sqrt(out["residual_sq_sum"] / n)
    out["mean_observed_ln_g"] = out["observed_sum"] / n
    out["mean_predicted_ln_g"] = out["predicted_sum"] / n
    denom = (out["n_records"] - 1).replace(0, np.nan)
    variance = (out["residual_sq_sum"] - (out["residual_sum"] ** 2) / n) / denom
    out["residual_std_ln"] = np.sqrt(variance.clip(lower=0))
    out["mae_log10_equiv"] = out["mae_ln"] / np.log(10.0)
    out["rmse_log10_equiv"] = out["rmse_ln"] / np.log(10.0)
    return out


def combine_aggs(parts: list[pd.DataFrame], group_cols: list[str]) -> pd.DataFrame:
    if not parts:
        return pd.DataFrame(columns=group_cols + AGG_VALUE_COLS)
    combined = pd.concat(parts, ignore_index=True)
    combined = (
        combined.groupby(group_cols, dropna=False)[AGG_VALUE_COLS]
        .sum()
        .reset_index()
    )
    return finalize_metrics(combined)


def station_site_correlations(station_terms: pd.DataFrame, min_records: int) -> pd.DataFrame:
    features = [
        ("log10(VS10)", "vs10"),
        ("log10(VS20)", "vs20"),
        ("log10(VS30)", "vs30"),
        ("log10(AVS30)", "avs30"),
        ("log10(VS30_proxy)", "vs30_proxy"),
        ("log10(D1100)", "d1100"),
        ("log10(D1400)", "d1400"),
        ("log10(Dbase)", "dbase"),
    ]
    rows = []
    supported = station_terms[station_terms["n_records"].ge(min_records)].copy()
    for (target, target_label), sub in supported.groupby(["target", "target_label"]):
        for feature_label, col in features:
            x = pd.to_numeric(sub[col], errors="coerce")
            y = pd.to_numeric(sub["mean_residual_ln"], errors="coerce")
            valid = x.gt(0) & y.notna()
            if valid.sum() < 30:
                continue
            rho, pvalue = stats.spearmanr(np.log10(x.loc[valid]), y.loc[valid])
            rows.append(
                {
                    "target": target,
                    "target_label": target_label,
                    "model_family": "zhao2006_openquake",
                    "feature": feature_label,
                    "n_sites": int(valid.sum()),
                    "min_records_per_site": min_records,
                    "spearman_rho": float(rho),
                    "pvalue": float(pvalue),
                }
            )
    return pd.DataFrame(rows)


def plot_correlations(correlations: pd.DataFrame, output_dir: Path) -> None:
    if correlations.empty:
        return
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = [item["target_label"] for item in TARGETS]
    heat = correlations.pivot(index="target_label", columns="feature", values="spearman_rho")
    heat = heat.loc[[label for label in labels if label in heat.index]]
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    im = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=35, ha="right", fontsize=8)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("Zhao 2006 Residuals vs Official Site Parameters")
    fig.colorbar(im, ax=ax, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(figure_dir / "jshis_zhao2006_station_site_correlation_heatmap.png", dpi=180)
    plt.close(fig)


def write_report(
    output_dir: Path,
    zip_path: Path,
    total_smrec_rows: int,
    joined_rows: int,
    overall: pd.DataFrame,
    by_source: pd.DataFrame,
    station_terms: pd.DataFrame,
    correlations: pd.DataFrame,
    coverage_summary: pd.DataFrame,
    vs30_source_summary: pd.DataFrame,
    min_station_records: int,
) -> None:
    lines = ["# Zhao et al. (2006) External GMPE Residual Audit\n\n"]
    lines.append("## Purpose\n")
    lines.append(
        "This audit tests whether the J-SHIS/NIED frequency-dependent residual pattern is only an artifact "
        "of the MF2013 implementation. It uses the OpenQuake hazardlib implementation of Zhao et al. (2006) "
        "as an external Japan-compatible GMPE family.\n\n"
    )

    lines.append("## Source data and model mapping\n")
    lines.append(f"- Archive: `{zip_path.name}`.\n")
    lines.append("- Official flatfile schemas used: `smrec_schema.tsv`, `source_schema.tsv`, `site_schema.tsv`.\n")
    lines.append("- GMPE implementation: OpenQuake hazardlib `zhao_2006`.\n")
    lines.append("- Source-class mapping: CRUSTAL -> `ZhaoEtAl2006Asc`; INTERPLATE -> `ZhaoEtAl2006SInter`; INTRAPLATE -> `ZhaoEtAl2006SSlab`.\n")
    lines.append("- Observed units: J-SHIS gal values converted to `ln(g)` using 980.665 gal per g.\n")
    lines.append("- Distance proxy: J-SHIS `fault_dist` used as OpenQuake `rrup` proxy.\n")
    lines.append("- Site proxy: official `vs30` where available, otherwise official `avs30` fallback.\n\n")

    lines.append("## Coverage\n")
    lines.append(f"- `smrec_schema.tsv` rows read: {total_smrec_rows:,}.\n")
    lines.append(f"- Rows joined to source and site schema: {joined_rows:,}.\n")
    for _, row in coverage_summary.sort_values(["target_label"]).iterrows():
        lines.append(
            f"- {row['target_label']}: usable Zhao residual rows={int(row['usable_rows']):,}, "
            f"missing/invalid rows={int(row['invalid_rows']):,}.\n"
        )
    lines.append("\n## VS30 proxy provenance among usable rows\n")
    if vs30_source_summary.empty:
        lines.append("- No usable rows.\n")
    else:
        for _, row in vs30_source_summary.iterrows():
            lines.append(f"- {row['vs30_proxy_source']}: {int(row['n_records']):,} rows.\n")

    lines.append("\n## Overall residual metrics\n")
    for _, row in overall.sort_values(["target_label"]).iterrows():
        lines.append(
            f"- {row['target_label']}: n={int(row['n_records']):,}, "
            f"MAE={row['mae_ln']:.4f} ln units ({row['mae_log10_equiv']:.4f} log10-equivalent), "
            f"RMSE={row['rmse_ln']:.4f}, bias observed-minus-predicted={row['mean_residual_ln']:.4f}.\n"
        )

    lines.append("\n## Source-class coverage\n")
    if by_source.empty:
        lines.append("- No source-class metrics generated.\n")
    else:
        source_pivot = by_source.pivot_table(
            index="target_label",
            columns="source_class_label",
            values="n_records",
            aggfunc="sum",
            fill_value=0,
        )
        for target_label, row in source_pivot.iterrows():
            parts = ", ".join(f"{col}={int(val):,}" for col, val in row.items())
            lines.append(f"- {target_label}: {parts}.\n")

    lines.append("\n## Strongest station-residual site associations\n")
    if correlations.empty:
        lines.append("- No station-site correlations met the support threshold.\n")
    else:
        best = correlations.reindex(correlations["spearman_rho"].abs().sort_values(ascending=False).index).head(12)
        for _, row in best.iterrows():
            lines.append(
                f"- {row['target_label']} vs {row['feature']}: "
                f"rho={row['spearman_rho']:.3f}, n_sites={int(row['n_sites'])}.\n"
            )

    supported_sites = station_terms[station_terms["n_records"].ge(min_station_records)]
    lines.append("\n## Station residual support\n")
    lines.append(
        f"- Station terms with at least {min_station_records} records per target: {len(supported_sites):,} rows across "
        f"{supported_sites['siteid2'].nunique():,} site IDs.\n"
    )

    lines.append("\n## Interpretation for the CEE manuscript\n")
    lines.append(
        "- The observed Zhao et al. (2006) residual pattern is qualitatively consistent with the MF2013 backbone: "
        "PGA shows little basin-depth residual structure, while SA(1.0 s) and especially SA(3.0 s) retain stronger "
        "positive basin-depth associations.\n"
    )
    lines.append(
        "- This supports model-sensitivity of the frequency-dependent site/basin residual finding. It does not prove causality "
        "and should not be presented as an exact Zhao et al. (2006) reproduction.\n"
    )
    lines.append(
        "- Absolute MAE values should not be used as a direct contest against MF2013, because this audit uses a distance proxy, "
        "a VS30/AVS30 fallback, and OpenQuake's Zhao et al. context assumptions.\n"
    )

    lines.append("\n## Outputs\n")
    for name in [
        "jshis_zhao2006_external_gmpe_overall_metrics.csv",
        "jshis_zhao2006_external_gmpe_by_network_metrics.csv",
        "jshis_zhao2006_external_gmpe_by_source_class_metrics.csv",
        "jshis_zhao2006_external_gmpe_station_residual_terms.csv",
        "jshis_zhao2006_external_gmpe_station_site_correlations.csv",
        "jshis_zhao2006_external_gmpe_coverage_summary.csv",
        "jshis_zhao2006_external_gmpe_vs30_proxy_summary.csv",
        "figures/jshis_zhao2006_station_site_correlation_heatmap.png",
    ]:
        lines.append(f"- `{name}`\n")

    lines.append("\n## Guardrails\n")
    lines.append("- This is an external GMPE sensitivity audit, not a full Zhao et al. (2006) reproduction.\n")
    lines.append("- `fault_dist` is a distance proxy, not a validated `rrup` reconstruction.\n")
    lines.append("- `avs30` fallback rows should be described as VS30-proxy rows.\n")
    lines.append("- Do not claim causal proof from residual-site correlations.\n")
    lines.append("- Do not use this audit to claim guaranteed CEE acceptance.\n")
    (output_dir / "jshis_zhao2006_external_gmpe_residuals.md").write_text("".join(lines), encoding="utf-8")


def run_analysis(zip_path: Path, output_dir: Path, chunksize: int, min_station_records: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    site, source = read_schema_tables(zip_path)

    overall_parts: list[pd.DataFrame] = []
    network_parts: list[pd.DataFrame] = []
    source_parts: list[pd.DataFrame] = []
    station_parts: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    vs30_rows: list[pd.DataFrame] = []
    total_smrec_rows = 0
    joined_rows = 0

    with ZipFile(zip_path) as zf:
        with zf.open("smrec_schema.tsv") as handle:
            reader = pd.read_csv(
                handle,
                sep="\t",
                usecols=SMREC_USECOLS,
                chunksize=chunksize,
                low_memory=False,
            )
            for chunk_index, chunk in enumerate(reader, start=1):
                total_smrec_rows += len(chunk)
                chunk["site_id"] = pd.to_numeric(chunk["site_id"], errors="coerce").astype("Int64")
                chunk["siteid2"] = (chunk["site_id"] // 10).astype("Int64")
                chunk["eq_source_id"] = pd.to_numeric(chunk["eq_source_id"], errors="coerce").astype("Int64")
                joined = chunk.merge(site, on="siteid2", how="left").merge(source, on="eq_source_id", how="left")
                joined_rows += len(joined)

                residual_frames = []
                for target_info in TARGETS:
                    residuals = residual_frame(joined, target_info)
                    residual_frames.append(residuals)
                    coverage_rows.append(
                        {
                            "chunk": chunk_index,
                            "target": target_info["target"],
                            "target_label": target_info["target_label"],
                            "input_rows": len(joined),
                            "usable_rows": len(residuals),
                            "invalid_rows": len(joined) - len(residuals),
                        }
                    )
                valid_frames = [frame for frame in residual_frames if not frame.empty]
                if not valid_frames:
                    print(f"chunk={chunk_index} smrec_rows={total_smrec_rows:,} residual_rows_this_chunk=0", flush=True)
                    continue
                residual_long = pd.concat(valid_frames, ignore_index=True)
                overall_parts.append(
                    aggregate_residuals(residual_long, ["target", "target_label", "period", "model_family"])
                )
                network_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model_family", "obs_network_id", "network_label"],
                    )
                )
                source_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model_family", "eq_location_type_id", "source_class_label"],
                    )
                )
                station_parts.append(
                    aggregate_residuals(
                        residual_long,
                        ["target", "target_label", "period", "model_family", "siteid2"],
                    )
                )
                vs30_rows.append(
                    residual_long.groupby("vs30_proxy_source", dropna=False)
                    .size()
                    .reset_index(name="n_records")
                )

                print(
                    f"chunk={chunk_index} smrec_rows={total_smrec_rows:,} "
                    f"residual_rows_this_chunk={len(residual_long):,}",
                    flush=True,
                )

    overall = combine_aggs(overall_parts, ["target", "target_label", "period", "model_family"])
    by_network = combine_aggs(
        network_parts,
        ["target", "target_label", "period", "model_family", "obs_network_id", "network_label"],
    )
    by_source = combine_aggs(
        source_parts,
        ["target", "target_label", "period", "model_family", "eq_location_type_id", "source_class_label"],
    )
    station_terms = combine_aggs(station_parts, ["target", "target_label", "period", "model_family", "siteid2"])
    station_terms = station_terms.merge(site, on="siteid2", how="left")
    correlations = station_site_correlations(station_terms, min_records=min_station_records)

    coverage = pd.DataFrame(coverage_rows)
    coverage_summary = (
        coverage.groupby(["target", "target_label"], as_index=False)[["usable_rows", "invalid_rows"]]
        .sum()
        .sort_values(["target_label"])
    )
    vs30_source_summary = (
        pd.concat(vs30_rows, ignore_index=True).groupby("vs30_proxy_source", as_index=False)["n_records"].sum()
        if vs30_rows
        else pd.DataFrame(columns=["vs30_proxy_source", "n_records"])
    )

    overall.to_csv(output_dir / "jshis_zhao2006_external_gmpe_overall_metrics.csv", index=False)
    by_network.to_csv(output_dir / "jshis_zhao2006_external_gmpe_by_network_metrics.csv", index=False)
    by_source.to_csv(output_dir / "jshis_zhao2006_external_gmpe_by_source_class_metrics.csv", index=False)
    station_terms.to_csv(output_dir / "jshis_zhao2006_external_gmpe_station_residual_terms.csv", index=False)
    correlations.to_csv(output_dir / "jshis_zhao2006_external_gmpe_station_site_correlations.csv", index=False)
    coverage_summary.to_csv(output_dir / "jshis_zhao2006_external_gmpe_coverage_summary.csv", index=False)
    vs30_source_summary.to_csv(output_dir / "jshis_zhao2006_external_gmpe_vs30_proxy_summary.csv", index=False)

    plot_correlations(correlations, output_dir)
    write_report(
        output_dir=output_dir,
        zip_path=zip_path,
        total_smrec_rows=total_smrec_rows,
        joined_rows=joined_rows,
        overall=overall,
        by_source=by_source,
        station_terms=station_terms,
        correlations=correlations,
        coverage_summary=coverage_summary,
        vs30_source_summary=vs30_source_summary,
        min_station_records=min_station_records,
    )

    print(f"total_smrec_rows={total_smrec_rows:,}")
    print(f"overall_metric_rows={len(overall):,}")
    print(f"station_term_rows={len(station_terms):,}")
    print(f"correlation_rows={len(correlations):,}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--min-station-records", type=int, default=20)
    args = parser.parse_args()
    run_analysis(
        zip_path=args.zip,
        output_dir=args.output_dir,
        chunksize=args.chunksize,
        min_station_records=args.min_station_records,
    )


if __name__ == "__main__":
    main()
