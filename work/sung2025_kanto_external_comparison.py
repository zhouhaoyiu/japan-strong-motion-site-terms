#!/usr/bin/env python3
"""Compare the Japanese SA(5 s) station field with Sung et al. (2025)."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
ARTICLE_DIR = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT_DIR = ARTICLE_DIR / "supplement"
FIGURE_DIR = ARTICLE_DIR / "figures"

DEFAULT_SUNG_CSV = (
    ROOT / "work" / "external_data" / "sung2025" / "BSSA-2024239_Supplement.csv"
)
FULL_PREDICTIONS = SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv"
SOURCE_PREDICTIONS = SUPPLEMENT_DIR / "jshis_source_conditioned_station_predictions.csv"

OUT_MATCHED = SUPPLEMENT_DIR / "sung2025_kanto_sa5_matched_station_fields.csv"
OUT_METRICS = SUPPLEMENT_DIR / "sung2025_kanto_sa5_external_metrics.csv"
OUT_BOOTSTRAP = SUPPLEMENT_DIR / "sung2025_kanto_sa5_spatial_block_bootstrap.csv"
OUT_METADATA = SUPPLEMENT_DIR / "sung2025_kanto_sa5_input_metadata.csv"
OUT_AUDIT = SUPPLEMENT_DIR / "sung2025_kanto_sa5_external_comparison.md"
OUT_FIGURE_PDF = FIGURE_DIR / "supplementary_figure_sung2025_kanto_external.pdf"
OUT_FIGURE_PNG = FIGURE_DIR / "supplementary_figure_sung2025_kanto_external.png"

PERIOD_S = 5.0
MATCH_SPACING_MULTIPLIER = 1.5
KANTO_CORRELATION_LENGTH_KM = 31.7
LN10 = float(np.log(10.0))

FIELD_LABELS = {
    "crustal_observed": "Crustal observed",
    "crustal_oof_spatial": "Crustal spatial OOF",
    "crustal_oof_physical": "Crustal physical OOF",
    "full_observed": "All-source observed",
    "full_oof_spatial": "All-source spatial OOF",
    "full_oof_physical": "All-source physical OOF",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def local_xy(
    lon: np.ndarray,
    lat: np.ndarray,
    lon0: float,
    lat0: float,
) -> tuple[np.ndarray, np.ndarray]:
    x = (lon - lon0) * 111.320 * np.cos(np.deg2rad(lat0))
    y = (lat - lat0) * 110.574
    return x, y


def load_sung_grid(path: Path) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    raw = pd.read_csv(path)
    expected = {"Lat", "Lon", "adjested_site"}
    if set(raw.columns) != expected:
        raise ValueError(f"unexpected Sung supplement columns: {list(raw.columns)}")
    grid = raw.rename(
        columns={
            "Lat": "sung_lat",
            "Lon": "sung_lon",
            "adjested_site": "sung_adjusted_site_ln",
        }
    ).copy()
    valid = grid[grid["sung_adjusted_site_ln"].notna()].copy().reset_index(names="grid_row")
    if len(grid) != 250_000 or valid.empty:
        raise ValueError("unexpected Sung supplement dimensions")

    lon0 = float(valid["sung_lon"].median())
    lat0 = float(valid["sung_lat"].median())
    valid["grid_x_km"], valid["grid_y_km"] = local_xy(
        valid["sung_lon"].to_numpy(float),
        valid["sung_lat"].to_numpy(float),
        lon0,
        lat0,
    )
    tree = cKDTree(valid[["grid_x_km", "grid_y_km"]].to_numpy(float))
    nearest = tree.query(valid[["grid_x_km", "grid_y_km"]].to_numpy(float), k=2)[0][
        :, 1
    ]
    spacing_km = float(np.median(nearest))
    metadata: dict[str, float | int | str] = {
        "input_path": str(path.relative_to(ROOT)),
        "input_bytes": path.stat().st_size,
        "input_sha256": sha256(path),
        "total_grid_rows": len(grid),
        "valid_grid_rows": len(valid),
        "projection_lon0": lon0,
        "projection_lat0": lat0,
        "median_valid_grid_spacing_km": spacing_km,
        "match_threshold_km": MATCH_SPACING_MULTIPLIER * spacing_km,
        "spatial_block_length_km": KANTO_CORRELATION_LENGTH_KM,
    }
    return valid, metadata


def extract_field(
    predictions: pd.DataFrame,
    field: str,
    value_column: str,
    model: str | None,
) -> pd.DataFrame:
    frame = predictions[predictions["period_s"].eq(PERIOD_S)].copy()
    if model is not None:
        frame = frame[frame["model"].eq(model)]
    columns = [
        "siteid2",
        "site_code",
        "network_label",
        "lon",
        "lat",
        "n_records",
        value_column,
    ]
    frame = frame[columns].dropna(subset=["lon", "lat", value_column]).copy()
    frame = frame.drop_duplicates("siteid2", keep="first")
    frame = frame.rename(columns={value_column: "station_value_log10"})
    frame["comparison_field"] = field
    frame["comparison_label"] = FIELD_LABELS[field]
    frame["station_value_ln"] = frame["station_value_log10"] * LN10
    return frame


def load_station_fields(
    full_path: Path,
    source_path: Path,
) -> pd.DataFrame:
    full = pd.read_csv(full_path)
    source = pd.read_csv(source_path)
    source = source[source["source_type"].eq("Crustal")].copy()
    blocks = [
        extract_field(full, "full_observed", "station_effect_log10", None),
        extract_field(
            full,
            "full_oof_spatial",
            "centered_oof_prediction_log10",
            "physical_spatial_hgb",
        ),
        extract_field(
            full,
            "full_oof_physical",
            "centered_oof_prediction_log10",
            "physical_hgb",
        ),
        extract_field(source, "crustal_observed", "station_effect_log10", None),
        extract_field(
            source,
            "crustal_oof_spatial",
            "centered_oof_prediction_log10",
            "physical_spatial_hgb",
        ),
        extract_field(
            source,
            "crustal_oof_physical",
            "centered_oof_prediction_log10",
            "physical_hgb",
        ),
    ]
    fields = pd.concat(blocks, ignore_index=True)
    if set(fields["comparison_field"]) != set(FIELD_LABELS):
        raise ValueError("one or more comparison fields are missing")
    return fields


def match_fields_to_grid(
    fields: pd.DataFrame,
    grid: pd.DataFrame,
    metadata: dict[str, float | int | str],
) -> pd.DataFrame:
    lon0 = float(metadata["projection_lon0"])
    lat0 = float(metadata["projection_lat0"])
    threshold = float(metadata["match_threshold_km"])
    tree = cKDTree(grid[["grid_x_km", "grid_y_km"]].to_numpy(float))
    matched_blocks = []
    for field, frame in fields.groupby("comparison_field", sort=False):
        frame = frame.copy()
        frame["station_x_km"], frame["station_y_km"] = local_xy(
            frame["lon"].to_numpy(float),
            frame["lat"].to_numpy(float),
            lon0,
            lat0,
        )
        distance, grid_index = tree.query(
            frame[["station_x_km", "station_y_km"]].to_numpy(float), k=1
        )
        frame["grid_match_distance_km"] = distance
        frame["matched_grid_index"] = grid_index
        frame = frame[frame["grid_match_distance_km"].le(threshold)].copy()
        grid_values = grid.iloc[frame["matched_grid_index"].to_numpy(int)]
        frame["sung_grid_row"] = grid_values["grid_row"].to_numpy(int)
        frame["sung_lon"] = grid_values["sung_lon"].to_numpy(float)
        frame["sung_lat"] = grid_values["sung_lat"].to_numpy(float)
        frame["sung_adjusted_site_ln"] = grid_values[
            "sung_adjusted_site_ln"
        ].to_numpy(float)
        frame["station_centered_ln"] = (
            frame["station_value_ln"] - frame["station_value_ln"].mean()
        )
        frame["sung_centered_ln"] = (
            frame["sung_adjusted_site_ln"] - frame["sung_adjusted_site_ln"].mean()
        )
        frame["spatial_block_x"] = np.floor(
            (frame["station_x_km"] - grid["grid_x_km"].min())
            / KANTO_CORRELATION_LENGTH_KM
        ).astype(int)
        frame["spatial_block_y"] = np.floor(
            (frame["station_y_km"] - grid["grid_y_km"].min())
            / KANTO_CORRELATION_LENGTH_KM
        ).astype(int)
        frame["spatial_block"] = (
            frame["spatial_block_x"].astype(str)
            + "_"
            + frame["spatial_block_y"].astype(str)
        )
        matched_blocks.append(frame)
    return pd.concat(matched_blocks, ignore_index=True)


def association_metrics(block: pd.DataFrame) -> dict[str, float | int]:
    observed = block["station_centered_ln"].to_numpy(float)
    external = block["sung_centered_ln"].to_numpy(float)
    pearson = float(np.corrcoef(observed, external)[0, 1])
    spearman = float(spearmanr(observed, external).statistic)
    slope = float(np.polyfit(external, observed, 1)[0])
    return {
        "n_stations": len(block),
        "n_unique_sung_cells": block["sung_grid_row"].nunique(),
        "n_spatial_blocks": block["spatial_block"].nunique(),
        "pearson": pearson,
        "spearman": spearman,
        "centered_rmse_ln": float(np.sqrt(np.mean((observed - external) ** 2))),
        "sign_agreement_fraction": float(np.mean(np.sign(observed) == np.sign(external))),
        "station_std_ln": float(np.std(observed, ddof=1)),
        "sung_std_ln": float(np.std(external, ddof=1)),
        "station_to_sung_std_ratio": float(
            np.std(observed, ddof=1) / np.std(external, ddof=1)
        ),
        "station_on_sung_slope": slope,
        "match_distance_km_q50": float(block["grid_match_distance_km"].quantile(0.5)),
        "match_distance_km_q95": float(block["grid_match_distance_km"].quantile(0.95)),
        "match_distance_km_max": float(block["grid_match_distance_km"].max()),
    }


def calculate_metrics(matched: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for field, block in matched.groupby("comparison_field", sort=False):
        row: dict[str, float | int | str | bool] = {
            "comparison_field": field,
            "comparison_label": FIELD_LABELS[field],
            "primary_source_matched_test": field == "crustal_observed",
        }
        row.update(association_metrics(block))
        rows.append(row)
    return pd.DataFrame(rows)


def spatial_block_bootstrap(
    matched: pd.DataFrame,
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    replicate_rows = []
    summary_rows = []
    for field_index, (field, block) in enumerate(
        matched.groupby("comparison_field", sort=False)
    ):
        block = block.reset_index(drop=True)
        groups = {
            key: indices.to_numpy(int)
            for key, indices in block.groupby("spatial_block", sort=True).groups.items()
        }
        group_keys = np.asarray(list(groups), dtype=object)
        rng = np.random.default_rng(seed + field_index)
        field_rows = []
        for replicate in range(replicates):
            sampled_keys = rng.choice(group_keys, size=len(group_keys), replace=True)
            sampled_index = np.concatenate([groups[key] for key in sampled_keys])
            sample = block.iloc[sampled_index].copy()
            sample["station_centered_ln"] = (
                sample["station_value_ln"] - sample["station_value_ln"].mean()
            )
            sample["sung_centered_ln"] = (
                sample["sung_adjusted_site_ln"]
                - sample["sung_adjusted_site_ln"].mean()
            )
            values = association_metrics(sample)
            field_rows.append(
                {
                    "comparison_field": field,
                    "replicate": replicate,
                    "sampled_spatial_blocks": len(group_keys),
                    "sampled_station_rows": len(sample),
                    "pearson": values["pearson"],
                    "spearman": values["spearman"],
                    "sign_agreement_fraction": values["sign_agreement_fraction"],
                    "station_on_sung_slope": values["station_on_sung_slope"],
                }
            )
        replicate_rows.extend(field_rows)
        field_frame = pd.DataFrame(field_rows)
        summary: dict[str, float | int | str] = {
            "comparison_field": field,
            "bootstrap_replicates": replicates,
        }
        for column in [
            "pearson",
            "spearman",
            "sign_agreement_fraction",
            "station_on_sung_slope",
        ]:
            summary[f"{column}_ci_low"] = float(field_frame[column].quantile(0.025))
            summary[f"{column}_ci_high"] = float(field_frame[column].quantile(0.975))
        summary_rows.append(summary)
    return pd.DataFrame(replicate_rows), pd.DataFrame(summary_rows)


def save_figure(
    grid: pd.DataFrame,
    matched: pd.DataFrame,
    metrics: pd.DataFrame,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.7,
            "savefig.dpi": 300,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.8), constrained_layout=True)
    global_sung = (
        grid["sung_adjusted_site_ln"] - grid["sung_adjusted_site_ln"].mean()
    )
    limit = float(np.quantile(np.abs(global_sung), 0.98))
    map_artist = axes[0, 0].scatter(
        grid["sung_lon"],
        grid["sung_lat"],
        c=global_sung,
        s=1.2,
        marker="s",
        linewidths=0,
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
        rasterized=True,
    )
    axes[0, 0].set(
        title="a  Sung et al. simulation-based field",
        xlabel="Longitude (degrees E)",
        ylabel="Latitude (degrees N)",
    )
    fig.colorbar(map_artist, ax=axes[0, 0], label="Centred adjusted site term (ln)")

    primary = matched[matched["comparison_field"].eq("crustal_observed")].copy()
    station_limit = float(np.quantile(np.abs(primary["station_centered_ln"]), 0.95))
    axes[0, 1].scatter(
        grid["sung_lon"],
        grid["sung_lat"],
        c="#D6D6D6",
        s=0.6,
        marker="s",
        linewidths=0,
        rasterized=True,
    )
    station_artist = axes[0, 1].scatter(
        primary["lon"],
        primary["lat"],
        c=primary["station_centered_ln"],
        s=17,
        edgecolors="white",
        linewidths=0.25,
        cmap="RdBu_r",
        vmin=-station_limit,
        vmax=station_limit,
        rasterized=True,
    )
    axes[0, 1].set(
        title="b  Source-matched empirical station field",
        xlabel="Longitude (degrees E)",
        ylabel="Latitude (degrees N)",
    )
    fig.colorbar(station_artist, ax=axes[0, 1], label="Centred crustal term (ln)")

    axes[1, 0].scatter(
        primary["sung_centered_ln"],
        primary["station_centered_ln"],
        s=13,
        alpha=0.28,
        color="#2B6CB0",
        edgecolors="none",
        rasterized=True,
    )
    primary["sung_bin"] = pd.qcut(
        primary["sung_centered_ln"], q=8, labels=False, duplicates="drop"
    )
    binned = primary.groupby("sung_bin").agg(
        x=("sung_centered_ln", "median"),
        y=("station_centered_ln", "median"),
        y25=("station_centered_ln", lambda values: values.quantile(0.25)),
        y75=("station_centered_ln", lambda values: values.quantile(0.75)),
    )
    axes[1, 0].errorbar(
        binned["x"],
        binned["y"],
        yerr=np.vstack([binned["y"] - binned["y25"], binned["y75"] - binned["y"]]),
        fmt="o",
        linestyle="none",
        markersize=4.2,
        capsize=2.5,
        color="#B22222",
        ecolor="#B22222",
        label="Equal-count-bin median and IQR",
        zorder=3,
    )
    primary_metric = metrics[metrics["comparison_field"].eq("crustal_observed")].iloc[0]
    axes[1, 0].text(
        0.03,
        0.97,
        f"Pearson r = {primary_metric['pearson']:.3f}\n"
        f"Spearman rho = {primary_metric['spearman']:.3f}\n"
        f"n = {int(primary_metric['n_stations'])}",
        transform=axes[1, 0].transAxes,
        va="top",
    )
    axes[1, 0].set(
        title="c  Source-matched station association",
        xlabel="Sung adjusted site term, centred (ln)",
        ylabel="Observed crustal station term, centred (ln)",
    )
    axes[1, 0].legend(frameon=False, fontsize=7.0, loc="lower right")

    shown = [
        "crustal_observed",
        "crustal_oof_spatial",
        "full_observed",
        "full_oof_spatial",
    ]
    plot_metrics = metrics.set_index("comparison_field").loc[shown].reset_index()
    positions = np.arange(len(plot_metrics))
    pearson_errors = np.vstack(
        [
            np.maximum(
                0.0, plot_metrics["pearson"] - plot_metrics["pearson_ci_low"]
            ),
            np.maximum(
                0.0, plot_metrics["pearson_ci_high"] - plot_metrics["pearson"]
            ),
        ]
    )
    spearman_errors = np.vstack(
        [
            np.maximum(
                0.0, plot_metrics["spearman"] - plot_metrics["spearman_ci_low"]
            ),
            np.maximum(
                0.0, plot_metrics["spearman_ci_high"] - plot_metrics["spearman"]
            ),
        ]
    )
    axes[1, 1].errorbar(
        positions - 0.07,
        plot_metrics["pearson"],
        yerr=pearson_errors,
        fmt="o",
        color="#2B6CB0",
        ecolor="#2B6CB0",
        capsize=3,
        label="Pearson",
    )
    axes[1, 1].errorbar(
        positions + 0.07,
        plot_metrics["spearman"],
        yerr=spearman_errors,
        fmt="s",
        color="#D55E00",
        ecolor="#D55E00",
        capsize=3,
        label="Spearman",
    )
    axes[1, 1].axhline(0.0, color="#555555", linestyle="--", linewidth=0.8)
    axes[1, 1].set_xticks(positions)
    axes[1, 1].set_xticklabels(
        [
            "Crustal\nobserved",
            "Crustal\nspatial OOF",
            "All-source\nobserved",
            "All-source\nspatial OOF",
        ]
    )
    axes[1, 1].set(
        title="d  External-field agreement",
        ylabel="Spatial association",
    )
    axes[1, 1].legend(frameon=False)
    for axis in axes[:1, :].flat:
        axis.set_aspect(1.0 / np.cos(np.deg2rad(float(grid["sung_lat"].median()))))
    for axis in axes.flat:
        axis.grid(alpha=0.15, linewidth=0.5)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audit(
    metrics: pd.DataFrame,
    metadata: dict[str, float | int | str],
) -> None:
    primary = metrics[metrics["comparison_field"].eq("crustal_observed")].iloc[0]
    lines = [
        "# Sung et al. (2025) Kanto SA(5 s) external comparison",
        "",
        "## Locked design",
        "",
        "- The external field is the unmodified valid portion of the publisher-supplied Kanto T=5 s adjusted-site-term CSV. No geographic response-value screen is applied.",
        f"- Each station is assigned its nearest valid simulation cell only when the distance is at most 1.5 times the median grid spacing ({float(metadata['match_threshold_km']):.3f} km).",
        "- Sung et al. report natural-log units. Japanese log10 station terms and out-of-fold predictions are multiplied by ln(10). Each matched field pair is centred independently because the two residual coordinates have arbitrary offsets.",
        "- The source-matched crustal observed field is the primary physical-comparability test. Crustal out-of-fold predictions, the all-source field and the physical-only model are retained as diagnostics.",
        f"- Spatial-block intervals use 2,000 replicates and the published Kanto VCM correlation length of {KANTO_CORRELATION_LENGTH_KM:g} km. Stations receive equal weight.",
        "- Matching and inference settings were fixed before calculating cross-field associations.",
        "",
        "## Results",
        "",
        f"The source-matched comparison contains {int(primary['n_stations']):,} stations in {int(primary['n_spatial_blocks'])} spatial blocks. Its Pearson correlation is {primary['pearson']:.3f} ({primary['pearson_ci_low']:.3f}--{primary['pearson_ci_high']:.3f}), and its Spearman correlation is {primary['spearman']:.3f} ({primary['spearman_ci_low']:.3f}--{primary['spearman_ci_high']:.3f}). Centred signs agree at {100.0 * primary['sign_agreement_fraction']:.1f}% of stations ({100.0 * primary['sign_agreement_fraction_ci_low']:.1f}%--{100.0 * primary['sign_agreement_fraction_ci_high']:.1f}%).",
        "",
        "| Field | Stations | Blocks | Pearson | 95% interval | Spearman | 95% interval | Sign agreement |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in metrics.iterrows():
        lines.append(
            f"| {row['comparison_label']} | {int(row['n_stations'])} | {int(row['n_spatial_blocks'])} | "
            f"{row['pearson']:.3f} | {row['pearson_ci_low']:.3f}--{row['pearson_ci_high']:.3f} | "
            f"{row['spearman']:.3f} | {row['spearman_ci_low']:.3f}--{row['spearman_ci_high']:.3f} | "
            f"{100.0 * row['sign_agreement_fraction']:.1f}% |"
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The publisher field combines the simulation-derived basin adjustment and spatial nonergodic site term. The empirical field is a station residual after the MF2013 AVS30 and D1400 terms and still contains average-path structure. Their association tests independent spatial concordance; it does not equate the two amplitudes or validate a production nonergodic GMM. The simulation uses crustal scenarios and one 3D velocity model, so the source-matched result is reported first and its spatial-block uncertainty is retained.",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_outputs(
    matched: pd.DataFrame,
    metrics: pd.DataFrame,
    metadata: dict[str, float | int | str],
) -> None:
    assert set(matched["comparison_field"]) == set(FIELD_LABELS)
    assert set(metrics["comparison_field"]) == set(FIELD_LABELS)
    assert int(metrics["primary_source_matched_test"].sum()) == 1
    assert matched["sung_adjusted_site_ln"].notna().all()
    assert matched["grid_match_distance_km"].max() <= float(
        metadata["match_threshold_km"]
    ) + 1e-12
    assert metrics[
        [
            "pearson",
            "spearman",
            "pearson_ci_low",
            "pearson_ci_high",
            "spearman_ci_low",
            "spearman_ci_high",
        ]
    ].notna().all().all()
    for path in [
        OUT_MATCHED,
        OUT_METRICS,
        OUT_BOOTSTRAP,
        OUT_METADATA,
        OUT_AUDIT,
        OUT_FIGURE_PDF,
        OUT_FIGURE_PNG,
    ]:
        assert path.exists() and path.stat().st_size > 0, path


def run(args: argparse.Namespace) -> None:
    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    grid, metadata = load_sung_grid(args.sung_csv)
    fields = load_station_fields(args.full_predictions, args.source_predictions)
    matched = match_fields_to_grid(fields, grid, metadata)
    metrics = calculate_metrics(matched)
    bootstrap_replicates, bootstrap_summary = spatial_block_bootstrap(
        matched, args.bootstrap_replicates, args.seed
    )
    metrics = metrics.merge(
        bootstrap_summary,
        on="comparison_field",
        how="left",
        validate="one_to_one",
    )
    metadata.update(
        {
            "source_doi": "https://doi.org/10.1785/0120240239",
            "source_archive_member": "BSSA-2024239_Supplement.csv",
            "published_units": "natural_log",
            "match_spacing_multiplier": MATCH_SPACING_MULTIPLIER,
            "bootstrap_replicates": args.bootstrap_replicates,
            "random_seed": args.seed,
        }
    )

    matched.to_csv(OUT_MATCHED, index=False)
    metrics.to_csv(OUT_METRICS, index=False)
    bootstrap_replicates.to_csv(OUT_BOOTSTRAP, index=False)
    pd.DataFrame([metadata]).to_csv(OUT_METADATA, index=False)
    save_figure(grid, matched, metrics)
    write_audit(metrics, metadata)
    validate_outputs(matched, metrics, metadata)
    primary = metrics[metrics["comparison_field"].eq("crustal_observed")].iloc[0]
    print(
        f"source-matched: n={int(primary['n_stations'])} "
        f"pearson={primary['pearson']:.3f} spearman={primary['spearman']:.3f}",
        flush=True,
    )
    print(f"wrote {OUT_AUDIT}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sung-csv", type=Path, default=DEFAULT_SUNG_CSV)
    parser.add_argument("--full-predictions", type=Path, default=FULL_PREDICTIONS)
    parser.add_argument("--source-predictions", type=Path, default=SOURCE_PREDICTIONS)
    parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=20260711)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
