#!/usr/bin/env python3
"""Build English Figure 9: representative spectra and uncertainty checks."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter


HERE = Path(__file__).resolve().parent
SUPP = HERE / "supplement"
FIGURES = HERE / "figures"

STATION_VALUES = SUPP / "jshis_multiperiod_official_response_station_values.csv"
BOOTSTRAP = SUPP / "jshis_official_hazard_response_bootstrap_ci.csv"
SURFACE_VALIDATION = SUPP / "jshis_multiperiod_station_surface_validation_summary.csv"
OUT_CASES = SUPP / "jshis_representative_city_spectrum_cases.csv"
OUT_UNCERTAINTY = SUPP / "jshis_representative_uncertainty_summary.csv"


CASES = [
    ("Tokyo", 35.6895, 139.6917),
    ("Osaka", 34.6937, 135.5023),
    ("Sapporo", 43.0618, 141.3545),
    ("Sendai", 38.2682, 140.8719),
    ("NE Hokkaido/Soya", 45.3000, 142.0000),
]


def haversine_km(lat1: float, lon1: float, lat2: pd.Series, lon2: pd.Series) -> pd.Series:
    radius = 6371.0088
    lat1_r = np.radians(lat1)
    lon1_r = np.radians(lon1)
    lat2_r = np.radians(lat2.astype(float))
    lon2_r = np.radians(lon2.astype(float))
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2.0) ** 2
    return pd.Series(radius * 2.0 * np.arcsin(np.sqrt(a)), index=lat2.index)


def nearest_sites(values: pd.DataFrame) -> pd.DataFrame:
    sites = values[["siteid2", "site_code", "network_label", "lon", "lat"]].drop_duplicates("siteid2").copy()
    rows: list[pd.Series] = []
    for case_name, case_lat, case_lon in CASES:
        dist = haversine_km(case_lat, case_lon, sites["lat"], sites["lon"])
        hit = sites.loc[dist.idxmin()].copy()
        hit["case_name"] = case_name
        hit["case_lat"] = case_lat
        hit["case_lon"] = case_lon
        hit["distance_to_case_km"] = float(dist.loc[hit.name])
        rows.append(hit)
    return pd.DataFrame(rows)


def build_case_table(values: pd.DataFrame) -> pd.DataFrame:
    picked = nearest_sites(values)
    q10 = values[values["probability_level"].eq("50y_10pct")].copy()
    out = picked.merge(q10, on=["siteid2", "site_code", "network_label", "lon", "lat"], how="left")
    out = out.sort_values(["case_name", "period_s"]).reset_index(drop=True)
    out["delta_pct_at_ordinate"] = 100.0 * (out["corrected_sa_g"] / out["official_sa_g"] - 1.0)
    keep = [
        "case_name",
        "case_lat",
        "case_lon",
        "site_code",
        "network_label",
        "lat",
        "lon",
        "distance_to_case_km",
        "period_s",
        "official_sa_g",
        "corrected_sa_g",
        "multiplier",
        "delta_pct_at_ordinate",
        "correction_source",
    ]
    return out[keep]


def build_uncertainty_table(bootstrap: pd.DataFrame, validation: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    b = bootstrap[bootstrap["probability_level"].eq("50y_10pct")].copy()
    for row in b.itertuples(index=False):
        rows.append(
            {
                "source": "station_bootstrap",
                "quantity": row.metric,
                "period_s": 3.0,
                "quantile": row.distribution_quantile,
                "observed": row.observed,
                "low": row.bootstrap_ci_low,
                "high": row.bootstrap_ci_high,
                "n": row.n_stations,
            }
        )

    v = validation[validation["method"].eq("site_space_hgb")].copy()
    for period_key, period_s in [("sa03", 0.3), ("sa10", 1.0), ("sa30", 3.0)]:
        sub = v[v["target_key"].eq(period_key)]
        vals = sub["rmse_reduction_vs_zero_pct"].astype(float)
        rows.append(
            {
                "source": "spatial_block_folds",
                "quantity": "rmse_reduction_vs_zero_pct",
                "period_s": period_s,
                "quantile": 0.50,
                "observed": float(vals.mean()),
                "low": float(vals.min()),
                "high": float(vals.max()),
                "n": int(sub["n_stations"].sum()),
            }
        )
    return pd.DataFrame(rows)


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "font.size": 10.4,
            "axes.titlesize": 11.4,
            "axes.labelsize": 10.2,
            "xtick.labelsize": 8.8,
            "ytick.labelsize": 8.8,
            "legend.fontsize": 8.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def panel(ax, label: str, title: str) -> None:
    ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=13, va="bottom")
    ax.text(-0.055, 1.04, title, transform=ax.transAxes, fontsize=11.2, va="bottom")


def fmt_period(value: float, _pos: int) -> str:
    return f"{value:g}" if value in {0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0} else ""


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    style()
    values = pd.read_csv(STATION_VALUES)
    cases = build_case_table(values)
    uncertainty = build_uncertainty_table(pd.read_csv(BOOTSTRAP), pd.read_csv(SURFACE_VALIDATION))
    cases.to_csv(OUT_CASES, index=False)
    uncertainty.to_csv(OUT_UNCERTAINTY, index=False)

    colors = {
        "Tokyo": "#4C78A8",
        "Osaka": "#F58518",
        "Sapporo": "#54A24B",
        "Sendai": "#B279A2",
        "NE Hokkaido/Soya": "#E45756",
    }

    fig = plt.figure(figsize=(12.8, 8.8))
    grid = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    for case_name, sub in cases.groupby("case_name", sort=False):
        sub = sub.sort_values("period_s")
        color = colors[case_name]
        label = f"{case_name} ({sub['site_code'].iloc[0]}, {sub['distance_to_case_km'].iloc[0]:.1f} km)"
        ax_a.plot(sub["period_s"], sub["official_sa_g"], color=color, lw=1.05, ls="--", alpha=0.72)
        ax_a.plot(sub["period_s"], sub["corrected_sa_g"], color=color, lw=1.8, marker="o", ms=3.0, label=label)
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5]))
    ax_a.xaxis.set_major_formatter(FuncFormatter(fmt_period))
    ax_a.set_xlabel("Period (s)")
    ax_a.set_ylabel("5%-damped SA, 50-year 10% (g)")
    ax_a.grid(color="#E6E6E6", which="both")
    panel(ax_a, "A", "Representative city-nearest spectra")
    ax_a.legend(frameon=False, loc="lower left")

    bar = cases[cases["period_s"].isin([1.0, 3.0, 5.0])].copy()
    labels = [case for case, *_ in CASES]
    x = np.arange(len(labels))
    width = 0.24
    for offset, period in [(-width, 1.0), (0.0, 3.0), (width, 5.0)]:
        vals = [
            float(bar[(bar["case_name"].eq(case)) & (bar["period_s"].eq(period))]["delta_pct_at_ordinate"].iloc[0])
            for case in labels
        ]
        ax_b.bar(x + offset, vals, width=width, label=f"{period:g} s")
    ax_b.axhline(0, color="#222222", lw=0.85)
    ax_b.set_xticks(x, labels, rotation=18, ha="right")
    ax_b.set_ylabel("Corrected minus official ordinate (%)")
    ax_b.grid(axis="y", color="#E6E6E6")
    panel(ax_b, "B", "Engineering-period ordinate changes")
    ax_b.legend(frameon=False, loc="best")

    b = uncertainty[
        uncertainty["source"].eq("station_bootstrap") & uncertainty["quantity"].eq("delta_pct")
    ].copy()
    rows = []
    for quantile, label in [(0.05, "5th percentile"), (0.50, "median"), (0.95, "95th percentile")]:
        row = b[b["quantile"].eq(quantile)].iloc[0]
        rows.append((label, row["observed"], row["low"], row["high"]))
    y = np.arange(len(rows))[::-1]
    for yi, (name, obs, lo, hi) in zip(y, rows, strict=True):
        ax_c.plot([lo, hi], [yi, yi], color="#4C78A8", lw=3.0, solid_capstyle="round")
        ax_c.scatter([obs], [yi], color="#222222", s=24, zorder=3)
    ax_c.set_yticks(y, [row[0] for row in rows])
    ax_c.axvline(0, color="#222222", lw=0.85)
    ax_c.set_xlabel("Bootstrap 95% interval for station-level SA(3.0 s) shift (%)")
    ax_c.grid(axis="x", color="#E6E6E6")
    panel(ax_c, "C", "Bootstrap uncertainty at SA(3.0 s)")

    v = uncertainty[uncertainty["source"].eq("spatial_block_folds")].copy().sort_values("period_s")
    ax_d.errorbar(
        v["period_s"].astype(str),
        v["observed"],
        yerr=[v["observed"] - v["low"], v["high"] - v["observed"]],
        fmt="o",
        color="#4C78A8",
        ecolor="#8FB5DA",
        capsize=4,
        lw=1.7,
    )
    ax_d.set_ylim(0, max(80, float(v["high"].max()) + 5))
    ax_d.set_xlabel("Period (s)")
    ax_d.set_ylabel("Spatial-block RMSE reduction (%)")
    ax_d.grid(axis="y", color="#E6E6E6")
    panel(ax_d, "D", "Fold range for correction-surface validation")

    fig.subplots_adjust(left=0.075, right=0.985, top=0.94, bottom=0.08)
    fig.savefig(FIGURES / "figure9_engineering_cases_uncertainty.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "figure9_engineering_cases_uncertainty.pdf", bbox_inches="tight")

    assert cases["case_name"].nunique() == len(CASES)
    assert cases.shape[0] == len(CASES) * 8
    assert uncertainty.shape[0] >= 15
    assert (FIGURES / "figure9_engineering_cases_uncertainty.pdf").exists()


if __name__ == "__main__":
    main()
