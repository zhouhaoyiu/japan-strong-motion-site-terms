#!/usr/bin/env python3
"""Propagate one frozen station field through official J-SHIS source-category maps."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

import jshis_event_adjusted_station_model as core


ROOT = Path(__file__).resolve().parents[1]
MAP_DIR = ROOT / "work" / "external_data" / "jshis_respmap"
PREDICTIONS = core.SUPPLEMENT_DIR / "jshis_event_adjusted_station_model_predictions.csv"

MAPS = {
    "all_earthquakes": (
        MAP_DIR / "P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip",
        "https://www.j-shis.bosai.go.jp/map/respmap/data/P/Y2020/"
        "P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip",
    ),
    "active_shallow": (
        MAP_DIR / "P-Y2020-RESP-MAP-AVR-LND_MTTL-T50.zip",
        "https://www.j-shis.bosai.go.jp/map/respmap/data/P/Y2020/"
        "P-Y2020-RESP-MAP-AVR-LND_MTTL-T50.zip",
    ),
    "subduction": (
        MAP_DIR / "P-Y2020-RESP-MAP-AVR-PPE_MTTL-T50.zip",
        "https://www.j-shis.bosai.go.jp/map/respmap/data/P/Y2020/"
        "P-Y2020-RESP-MAP-AVR-PPE_MTTL-T50.zip",
    ),
}

OUT_VALUES = core.SUPPLEMENT_DIR / "jshis_source_category_surface_spectrum_values.csv"
OUT_SUMMARY = core.SUPPLEMENT_DIR / "jshis_source_category_surface_spectrum_summary.csv"
OUT_STATIONS = core.SUPPLEMENT_DIR / "jshis_source_category_sa3_station_comparison.csv"
OUT_AUDIT = core.SUPPLEMENT_DIR / "jshis_source_category_response_sensitivity.md"
OUT_FIGURE_PDF = core.FIGURE_DIR / "supplementary_figure_source_category_response.pdf"
OUT_FIGURE_PNG = core.FIGURE_DIR / "supplementary_figure_source_category_response.png"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_values() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    predictions = pd.read_csv(PREDICTIONS, dtype={"meshcode250_str": "string"})
    coefficients = core.load_coefficients(core.DEFAULT_COEFFICIENTS)
    value_blocks = []
    summary_blocks = []
    hashes = {}
    for category, (path, _) in MAPS.items():
        if not path.exists():
            raise FileNotFoundError(path)
        hashes[category] = sha256(path)
        values, summary = core.build_hazard_values(predictions, coefficients, path)
        values["source_category"] = category
        summary["source_category"] = category
        value_blocks.append(
            values[
                [
                    "source_category", "siteid2", "site_code", "network_label",
                    "meshcode250_str", "lon", "lat", "period_s", "period_code",
                    "probability_level", "official_vs400_sa_g", "surface_reference_factor",
                    "ergodic_surface_sa_g", "oof_multiplier", "adjusted_surface_sa_g",
                ]
            ]
        )
        summary_blocks.append(summary)
        print(f"source category: {category} rows={len(values):,}", flush=True)
    return (
        pd.concat(value_blocks, ignore_index=True),
        pd.concat(summary_blocks, ignore_index=True),
        hashes,
    )


def station_comparison(values: pd.DataFrame) -> pd.DataFrame:
    selected = values[
        values["period_s"].eq(3.0) & values["probability_level"].eq("50y_10pct")
    ].copy()
    identity = [
        "siteid2", "site_code", "network_label", "meshcode250_str", "lon", "lat",
        "surface_reference_factor", "oof_multiplier",
    ]
    wide = selected.pivot(
        index=identity,
        columns="source_category",
        values=["official_vs400_sa_g", "ergodic_surface_sa_g", "adjusted_surface_sa_g"],
    )
    wide.columns = [f"{measure}_{category}" for measure, category in wide.columns]
    wide = wide.reset_index()
    shallow = wide["official_vs400_sa_g_active_shallow"].to_numpy(float)
    subduction = wide["official_vs400_sa_g_subduction"].to_numpy(float)
    wide["dominant_source_category"] = np.where(
        subduction >= shallow, "subduction", "active_shallow"
    )
    floor = np.finfo(float).tiny
    wide["log10_subduction_to_shallow_ratio"] = np.log10(
        np.maximum(subduction, floor) / np.maximum(shallow, floor)
    )
    return wide


def plot_results(values: pd.DataFrame, stations: pd.DataFrame) -> None:
    selected = values[
        values["period_s"].eq(3.0) & values["probability_level"].eq("50y_10pct")
    ]
    colors = {"active_shallow": "#D55E00", "subduction": "#0072B2"}
    labels = {"active_shallow": "Active shallow", "subduction": "Subduction"}
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), constrained_layout=True)

    ax = axes[0, 0]
    points = ax.scatter(
        stations["lon"], stations["lat"],
        c=stations["log10_subduction_to_shallow_ratio"], cmap="RdBu_r",
        vmin=-1.5, vmax=1.5, s=13, linewidths=0, alpha=0.8,
    )
    ax.set(
        xlabel="Longitude (degrees E)", ylabel="Latitude (degrees N)",
        title="a  Source-category contrast at SA(3.0 s)",
    )
    fig.colorbar(points, ax=ax, label="log$_{10}$(subduction / shallow)")

    ax = axes[0, 1]
    x = stations["official_vs400_sa_g_active_shallow"].to_numpy(float)
    y = stations["official_vs400_sa_g_subduction"].to_numpy(float)
    positive = (x > 0) & (y > 0)
    ax.scatter(x[positive], y[positive], s=12, alpha=0.35, color="#555555", linewidths=0)
    limits = [min(x[positive].min(), y[positive].min()), max(x.max(), y.max())]
    ax.plot(limits, limits, "--", color="#222222", lw=1)
    ax.text(0.04, 0.96, "Subduction larger", transform=ax.transAxes, va="top", color="#0072B2")
    ax.text(0.96, 0.04, "Active shallow larger", transform=ax.transAxes, ha="right", color="#D55E00")
    ax.set(
        xscale="log", yscale="log", xlim=limits, ylim=limits,
        xlabel="Active-shallow SA(3.0 s) (g)", ylabel="Subduction SA(3.0 s) (g)",
        title="b  Source-category contrast at matched stations",
    )

    ax = axes[1, 0]
    groups = [
        stations.loc[stations["dominant_source_category"].eq(category), "oof_multiplier"]
        for category in ["active_shallow", "subduction"]
    ]
    boxes = ax.boxplot(groups, tick_labels=["Active shallow", "Subduction"], patch_artist=True)
    for box, color in zip(boxes["boxes"], [colors["active_shallow"], colors["subduction"]], strict=True):
        box.set(facecolor=color, alpha=0.55)
    ax.axhline(1.0, color="#333333", ls="--", lw=1)
    ax.set(ylabel="Cross-validated station multiplier", title="c  Multiplier by dominant category")

    ax = axes[1, 1]
    probability = values[values["probability_level"].eq("50y_10pct")]
    for category in ["active_shallow", "subduction"]:
        block = probability[probability["source_category"].eq(category)]
        median = block.groupby("period_s", as_index=False)[
            ["ergodic_surface_sa_g", "adjusted_surface_sa_g"]
        ].median()
        ax.plot(
            median["period_s"], median["ergodic_surface_sa_g"], "o--",
            color=colors[category], alpha=0.65,
        )
        ax.plot(
            median["period_s"], median["adjusted_surface_sa_g"], "o-",
            color=colors[category], label=labels[category],
        )
    ax.set(
        xscale="log", yscale="log", xlabel="Period (s)", ylabel="Median surface SA (g)",
        title="d  Source-category surface spectra",
    )
    source_legend = ax.legend(frameon=False, loc="upper right")
    ax.add_artist(source_legend)
    ax.legend(
        handles=[
            Line2D([0], [0], color="#666666", ls="--", marker="o", label="Surface reference"),
            Line2D([0], [0], color="#666666", ls="-", marker="o", label="Station adjusted"),
        ],
        frameon=False,
        loc="lower left",
        fontsize=8,
    )
    for axis in axes.flat:
        axis.grid(alpha=0.2, linewidth=0.6)
    fig.savefig(OUT_FIGURE_PDF, bbox_inches="tight")
    fig.savefig(OUT_FIGURE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_audit(values: pd.DataFrame, stations: pd.DataFrame, hashes: dict[str, str]) -> None:
    dominant = stations["dominant_source_category"].value_counts()
    sa3 = values[
        values["period_s"].eq(3.0)
        & values["probability_level"].eq("50y_10pct")
        & values["source_category"].isin(["active_shallow", "subduction"])
    ]
    lines = [
        "# J-SHIS source-category response-spectrum sensitivity",
        "",
        "## Data boundary",
        "",
        "- The three inputs are official J-SHIS 2020 engineering-bedrock response-spectrum maps for all earthquakes, active shallow earthquakes and subduction earthquakes.",
    ]
    for category, (path, url) in MAPS.items():
        lines.append(f"- {category}: `{path.name}`; SHA-256 `{hashes[category]}`; {url}")
    lines.extend(
        [
            "- Source-category ordinates at a fixed exceedance probability are evaluated separately. They are not summed, because probability aggregation is nonlinear.",
            "- The same frozen spatial-block station multiplier is applied to each product. This isolates sensitivity to the official source category; it does not estimate source-specific station terms.",
            "",
            "## SA(3.0 s), 10% in 50 years",
            "",
            f"- Paired stations: {len(stations):,}.",
            f"- Active-shallow ordinate is larger at {int(dominant.get('active_shallow', 0)):,} stations; subduction ordinate is at least as large at {int(dominant.get('subduction', 0)):,} stations.",
        ]
    )
    for category, block in sa3.groupby("source_category", sort=True):
        lines.append(
            f"- {category}: median official ordinate {block['official_vs400_sa_g'].median():.3f} g; "
            f"median surface reference {block['ergodic_surface_sa_g'].median():.3f} g; "
            f"median station-adjusted value {block['adjusted_surface_sa_g'].median():.3f} g."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The calculation demonstrates the station correction on official maps spanning distinct source regimes. It remains a matched-station response-spectrum sensitivity. A production J-SHIS PSHA rerun would require the complete official source occurrence model, logic-tree weights and implementation details inside the hazard integral.",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate(values: pd.DataFrame, summary: pd.DataFrame, stations: pd.DataFrame) -> None:
    assert set(values["source_category"]) == set(MAPS)
    assert set(summary["source_category"]) == set(MAPS)
    assert set(values["period_s"]) == {spec.period_s for spec in core.PERIODS}
    assert len(stations) > 1_000
    assert stations[["official_vs400_sa_g_active_shallow", "official_vs400_sa_g_subduction"]].notna().all().all()
    for path in [OUT_VALUES, OUT_SUMMARY, OUT_STATIONS, OUT_AUDIT, OUT_FIGURE_PDF, OUT_FIGURE_PNG]:
        assert path.exists() and path.stat().st_size > 0, path


def run() -> None:
    values, summary, hashes = build_values()
    stations = station_comparison(values)
    values.to_csv(OUT_VALUES, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    stations.to_csv(OUT_STATIONS, index=False)
    plot_results(values, stations)
    write_audit(values, stations, hashes)
    validate(values, summary, stations)
    print(f"wrote {OUT_AUDIT}", flush=True)


if __name__ == "__main__":
    sys.exit(run())
