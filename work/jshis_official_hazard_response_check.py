#!/usr/bin/env python3
"""Check SA(3.0 s) station-correction impacts against official J-SHIS products.

The station model provides a log10 SA(3.0 s) correction at held-out stations.
This script joins that correction to public J-SHIS site mesh codes and applies
it to the official 2020 response-spectrum map ordinates at the same 250 m
mesh. It also stores official J-SHIS hazard-curve API outputs for example
stations as a reference for local hazard-curve shape.

The calculation changes only the SA(3.0 s) ordinate. It does not replace a full
PSHA calculation with source recurrence, logic-tree, or model-uncertainty
propagation.
"""

from __future__ import annotations

import argparse
import json
import math
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedLocator, FuncFormatter


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IMPACT = PROJECT_ROOT / "outputs" / "jshis_station_hazard_impact_by_station.csv"
DEFAULT_EXAMPLES = PROJECT_ROOT / "outputs" / "jshis_station_hazard_impact_examples.csv"
DEFAULT_SITE_SCHEMA = PROJECT_ROOT / "outputs" / "jshis_site_schema_v2024_sub1.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures"
DEFAULT_RESP_ZIP = PROJECT_ROOT / "work" / "external_data" / "jshis_respmap" / "P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip"
DEFAULT_API_CACHE = PROJECT_ROOT / "work" / "external_data" / "jshis_hazard_api_cache"

G_IN_CM_S2 = 980.665
PROBABILITY_COLUMNS = {
    "T50_P02_SA": "50y_2pct",
    "T50_P05_SA": "50y_5pct",
    "T50_P10_SA": "50y_10pct",
    "T50_P39_SA": "50y_39pct",
}
PERIOD_CODES = {
    "P0010": 0.1,
    "P0020": 0.2,
    "P0030": 0.3,
    "P0050": 0.5,
    "P0100": 1.0,
    "P0200": 2.0,
    "P0300": 3.0,
    "P0500": 5.0,
}
EXAMPLE_ORDER = ["maximum_negative", "q05", "q25", "q50", "q75", "q95", "maximum_positive"]
API_EXAMPLE_TYPES = ["maximum_negative", "q50", "maximum_positive"]


def normalize_mesh(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def load_station_frame(impact_path: Path, site_schema_path: Path) -> pd.DataFrame:
    impact = pd.read_csv(impact_path)
    schema = pd.read_csv(
        site_schema_path,
        usecols=["siteid2", "site_code", "lon", "lat", "meshcode250", "meshcode3"],
    )
    frame = impact.merge(schema, on=["siteid2", "site_code"], how="left", validate="many_to_one")
    frame["meshcode250_str"] = frame["meshcode250"].map(normalize_mesh)
    frame = frame.dropna(subset=["lon", "lat"])
    frame = frame[frame["meshcode250_str"].ne("")]
    return frame.reset_index(drop=True)


def match_examples(examples_path: Path, station_frame: pd.DataFrame) -> pd.DataFrame:
    examples = pd.read_csv(examples_path)
    rows = []
    for _, example in examples.iterrows():
        candidates = station_frame[station_frame["site_code"].eq(example["site_code"])].copy()
        if candidates.empty:
            continue
        diff = (candidates["sa3_uhs_multiplier"].astype(float) - float(example["sa3_uhs_multiplier"])).abs()
        row = candidates.iloc[int(diff.to_numpy().argmin())].copy()
        row["example_type"] = example["example_type"]
        rows.append(row)
    out = pd.DataFrame(rows)
    out["example_type"] = pd.Categorical(out["example_type"], EXAMPLE_ORDER, ordered=True)
    return out.sort_values("example_type").reset_index(drop=True)


def fetch_jshis_hazard_curve(
    row: pd.Series,
    period_years: int,
    api_cache: Path,
    timeout: int = 30,
    sleep_seconds: float = 0.2,
) -> pd.DataFrame:
    api_cache.mkdir(parents=True, exist_ok=True)
    site_code = str(row["site_code"])
    cache_file = api_cache / f"{site_code}_{int(row['siteid2'])}_Y2024_AVR_TTL_MTTL_T{period_years}.json"
    if cache_file.exists():
        data = json.loads(cache_file.read_text(encoding="utf-8"))
    else:
        query = urllib.parse.urlencode({"position": f"{row['lon']:.6f},{row['lat']:.6f}", "epsg": "4326"})
        url = f"https://www.j-shis.bosai.go.jp/map/api/pshm/Y2024/AVR/TTL_MTTL/T{period_years}/hzcv.json?{query}"
        request = urllib.request.Request(url, headers={"User-Agent": "station-residual-reproducibility-check/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = json.loads(response.read().decode("utf-8"))
        cache_file.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
        time.sleep(sleep_seconds)

    status = data.get("status", "")
    if status != "Success":
        raise RuntimeError(f"J-SHIS API returned {status} for {site_code} T{period_years}")

    values = [float(v) for v in data["sim"]["value"]]
    probabilities = [float(v) for v in data["prob"]["value"]]
    rows = []
    for sim_value, probability in zip(values, probabilities, strict=True):
        if probability >= 1.0:
            annual_rate = math.inf
        elif probability <= 0.0:
            annual_rate = 0.0
        else:
            annual_rate = -math.log1p(-probability) / float(period_years)
        rows.append(
            {
                "example_type": row["example_type"],
                "siteid2": int(row["siteid2"]),
                "site_code": site_code,
                "lon": float(row["lon"]),
                "lat": float(row["lat"]),
                "meshcode3": normalize_mesh(row["meshcode3"]),
                "meshcode250": row["meshcode250_str"],
                "period_years": period_years,
                "sim_type": data["sim"].get("type", ""),
                "sim_unit": data["sim"].get("unit", ""),
                "intensity_value": sim_value,
                "exceedance_probability": probability,
                "annual_rate": annual_rate,
            }
        )
    return pd.DataFrame(rows)


def build_api_reference(examples: pd.DataFrame, api_cache: Path) -> pd.DataFrame:
    frames = []
    api_examples = examples[examples["example_type"].astype(str).isin(API_EXAMPLE_TYPES)].copy()
    for _, row in api_examples.iterrows():
        for period_years in [30, 50]:
            frames.append(fetch_jshis_hazard_curve(row, period_years, api_cache))
    return pd.concat(frames, ignore_index=True)


def find_period_entry(zf: zipfile.ZipFile, period_code: str) -> str | None:
    code = period_code.upper()
    candidates = []
    for name in zf.namelist():
        upper = Path(name).name.upper()
        if not upper.endswith(".CSV"):
            continue
        if code not in upper:
            continue
        candidates.append(name)
    if not candidates:
        return None
    candidates = sorted(candidates, key=lambda item: (("PSV" in item.upper()), len(item), item))
    return candidates[0]


def read_response_period(zf: zipfile.ZipFile, period_code: str, meshcodes: set[str]) -> pd.DataFrame:
    entry = find_period_entry(zf, period_code)
    if entry is None:
        raise FileNotFoundError(f"No response-spectrum CSV found for {period_code}")
    names = ["CODE", "T50_P02_SA", "T50_P05_SA", "T50_P10_SA", "T50_P39_SA"]
    selected = []
    with zf.open(entry) as handle:
        chunks = pd.read_csv(handle, comment="#", header=None, names=names, dtype=str, chunksize=500_000)
        for chunk in chunks:
            chunk = chunk[chunk["CODE"].notna()].copy()
            chunk = chunk[~chunk["CODE"].astype(str).str.upper().eq("CODE")]
            chunk["meshcode250_str"] = chunk["CODE"].map(normalize_mesh)
            chunk = chunk[chunk["meshcode250_str"].isin(meshcodes)]
            if chunk.empty:
                continue
            for col in names[1:]:
                chunk[col] = pd.to_numeric(chunk[col], errors="coerce")
            chunk["period_code"] = period_code
            chunk["period_s"] = PERIOD_CODES[period_code]
            selected.append(chunk[["meshcode250_str", "period_code", "period_s", *names[1:]]])
    if not selected:
        return pd.DataFrame(columns=["meshcode250_str", "period_code", "period_s", *names[1:]])
    return pd.concat(selected, ignore_index=True)


def load_response_maps(resp_zip: Path, meshcodes: set[str]) -> pd.DataFrame:
    if not resp_zip.exists():
        raise FileNotFoundError(resp_zip)
    with zipfile.ZipFile(resp_zip) as zf:
        frames = [read_response_period(zf, code, meshcodes) for code in PERIOD_CODES]
    return pd.concat(frames, ignore_index=True)


def build_sa3_station_values(station_frame: pd.DataFrame, response_maps: pd.DataFrame) -> pd.DataFrame:
    sa3_map = response_maps[response_maps["period_code"].eq("P0300")].copy()
    merged = station_frame.merge(sa3_map, on="meshcode250_str", how="left", validate="many_to_one")
    rows = []
    for _, row in merged.iterrows():
        for source_col, label in PROBABILITY_COLUMNS.items():
            official_cm_s2 = row.get(source_col)
            if pd.isna(official_cm_s2):
                continue
            official_g = float(official_cm_s2) / G_IN_CM_S2
            factor = float(row["sa3_uhs_multiplier"])
            rows.append(
                {
                    "siteid2": int(row["siteid2"]),
                    "site_code": row["site_code"],
                    "network_label": row["network_label"],
                    "lon": float(row["lon"]),
                    "lat": float(row["lat"]),
                    "meshcode250": row["meshcode250_str"],
                    "probability_level": label,
                    "official_sa3_cm_s2": float(official_cm_s2),
                    "official_sa3_g": official_g,
                    "station_multiplier": factor,
                    "corrected_sa3_g": official_g * factor,
                    "delta_sa3_g": official_g * (factor - 1.0),
                    "delta_pct": 100.0 * (factor - 1.0),
                    "sa3_correction_log10": float(row["sa3_correction_log10"]),
                    "n_records": int(row["n_records"]),
                }
            )
    return pd.DataFrame(rows)


def build_uhs_examples(examples: pd.DataFrame, response_maps: pd.DataFrame) -> pd.DataFrame:
    long_map = response_maps.melt(
        id_vars=["meshcode250_str", "period_code", "period_s"],
        value_vars=list(PROBABILITY_COLUMNS),
        var_name="probability_column",
        value_name="official_sa_cm_s2",
    )
    long_map["probability_level"] = long_map["probability_column"].map(PROBABILITY_COLUMNS)
    merged = examples.merge(long_map, on="meshcode250_str", how="left", validate="many_to_many")
    merged = merged.dropna(subset=["official_sa_cm_s2"]).copy()
    merged["official_sa_g"] = merged["official_sa_cm_s2"].astype(float) / G_IN_CM_S2
    merged["station_corrected_sa_g"] = merged["official_sa_g"]
    is_sa3 = np.isclose(merged["period_s"].astype(float), 3.0)
    merged.loc[is_sa3, "station_corrected_sa_g"] = (
        merged.loc[is_sa3, "official_sa_g"].astype(float) * merged.loc[is_sa3, "sa3_uhs_multiplier"].astype(float)
    )
    return merged[
        [
            "example_type",
            "siteid2",
            "site_code",
            "network_label",
            "meshcode250_str",
            "period_code",
            "period_s",
            "probability_level",
            "official_sa_g",
            "station_corrected_sa_g",
            "sa3_uhs_multiplier",
            "sa3_correction_log10",
        ]
    ].rename(columns={"meshcode250_str": "meshcode250"})


def build_summary(sa3_values: pd.DataFrame, examples: pd.DataFrame, api_reference: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for level, sub in sa3_values.groupby("probability_level", sort=False):
        rows.append(
            {
                "metric": f"official_sa3_{level}",
                "n_station_values": len(sub),
                "official_sa3_g_q05": sub["official_sa3_g"].quantile(0.05),
                "official_sa3_g_q50": sub["official_sa3_g"].quantile(0.50),
                "official_sa3_g_q95": sub["official_sa3_g"].quantile(0.95),
                "corrected_sa3_g_q05": sub["corrected_sa3_g"].quantile(0.05),
                "corrected_sa3_g_q50": sub["corrected_sa3_g"].quantile(0.50),
                "corrected_sa3_g_q95": sub["corrected_sa3_g"].quantile(0.95),
                "delta_pct_q05": sub["delta_pct"].quantile(0.05),
                "delta_pct_q50": sub["delta_pct"].quantile(0.50),
                "delta_pct_q95": sub["delta_pct"].quantile(0.95),
            }
        )
    multiplier = examples["sa3_uhs_multiplier"].astype(float)
    rows.append(
        {
            "metric": "example_station_multipliers",
            "n_station_values": len(examples),
            "official_sa3_g_q05": np.nan,
            "official_sa3_g_q50": np.nan,
            "official_sa3_g_q95": np.nan,
            "corrected_sa3_g_q05": np.nan,
            "corrected_sa3_g_q50": np.nan,
            "corrected_sa3_g_q95": np.nan,
            "delta_pct_q05": 100.0 * (multiplier.quantile(0.05) - 1.0),
            "delta_pct_q50": 100.0 * (multiplier.quantile(0.50) - 1.0),
            "delta_pct_q95": 100.0 * (multiplier.quantile(0.95) - 1.0),
        }
    )
    finite_api = api_reference[np.isfinite(api_reference["annual_rate"].astype(float))].copy()
    rows.append(
        {
            "metric": "official_api_curve_points",
            "n_station_values": len(finite_api),
            "official_sa3_g_q05": np.nan,
            "official_sa3_g_q50": np.nan,
            "official_sa3_g_q95": np.nan,
            "corrected_sa3_g_q05": np.nan,
            "corrected_sa3_g_q50": np.nan,
            "corrected_sa3_g_q95": np.nan,
            "delta_pct_q05": np.nan,
            "delta_pct_q50": np.nan,
            "delta_pct_q95": np.nan,
        }
    )
    return pd.DataFrame(rows)


def write_report(
    summary: pd.DataFrame,
    examples: pd.DataFrame,
    sa3_values: pd.DataFrame,
    uhs_examples: pd.DataFrame,
    output_path: Path,
) -> None:
    def markdown_table(frame: pd.DataFrame) -> str:
        columns = list(frame.columns)
        lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
        for _, item in frame.iterrows():
            vals = []
            for col in columns:
                value = item[col]
                if isinstance(value, (float, np.floating)):
                    vals.append("" if pd.isna(value) else f"{float(value):.3f}")
                else:
                    vals.append("" if pd.isna(value) else str(value))
            lines.append("| " + " | ".join(vals) + " |")
        return "\n".join(lines)

    median = sa3_values[sa3_values["probability_level"].eq("50y_10pct")]
    example_sa3 = uhs_examples[
        uhs_examples["probability_level"].eq("50y_10pct") & np.isclose(uhs_examples["period_s"].astype(float), 3.0)
    ].copy()
    lines = [
        "# Official J-SHIS Hazard and Response-Spectrum Check",
        "",
        "This check applies the held-out SA(3.0 s) station correction to official J-SHIS 2020 response-spectrum map ordinates at the same 250 m station mesh. The response-spectrum data are on engineering bedrock (Vs=400 m/s) and use 5% damping. The calculation adjusts only the SA(3.0 s) ordinate and does not constitute a full PSHA rerun.",
        "",
        "## Data sources",
        "",
        "- J-SHIS PSHM hazard-curve API, Y2024, average case, all earthquakes, 30-year and 50-year periods.",
        "- J-SHIS response-spectrum map data, Y2020, average case, all earthquakes, 50-year period.",
        "",
        "## Station-level official SA(3.0 s) map values",
        "",
        f"- Matched station-probability rows: {len(sa3_values)}",
        f"- Matched held-out stations at 50-year 10% level: {len(median)}",
        f"- 50-year 10% official SA(3.0 s) q05/q50/q95: {median['official_sa3_g'].quantile(0.05):.3f} / {median['official_sa3_g'].quantile(0.50):.3f} / {median['official_sa3_g'].quantile(0.95):.3f} g",
        f"- 50-year 10% corrected SA(3.0 s) q05/q50/q95: {median['corrected_sa3_g'].quantile(0.05):.3f} / {median['corrected_sa3_g'].quantile(0.50):.3f} / {median['corrected_sa3_g'].quantile(0.95):.3f} g",
        "",
        "## Example stations at 50-year 10%",
        "",
        "| example_type | site_code | meshcode250 | official_SA3_g | corrected_SA3_g | multiplier |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in example_sa3.sort_values("example_type").iterrows():
        lines.append(
            f"| {row['example_type']} | {row['site_code']} | {row['meshcode250']} | "
            f"{row['official_sa_g']:.3f} | {row['station_corrected_sa_g']:.3f} | {row['sa3_uhs_multiplier']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Summary table",
            "",
            markdown_table(summary),
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def build_figure(
    station_frame: pd.DataFrame,
    examples: pd.DataFrame,
    api_reference: pd.DataFrame,
    sa3_values: pd.DataFrame,
    uhs_examples: pd.DataFrame,
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig = plt.figure(figsize=(12.8, 8.6))
    grid = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    factors = station_frame["sa3_uhs_multiplier"].astype(float)
    ax_a.hist(factors, bins=np.geomspace(max(factors.min(), 0.08), factors.max() * 1.05, 34), color="#4C78A8", edgecolor="white")
    ax_a.set_xscale("log")
    ax_a.set_xlim(0.1, 1.8)
    ax_a.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1.0, 1.5]))
    ax_a.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax_a.axvline(1.0, color="black", lw=1.0)
    for q, ls in [(0.05, "--"), (0.50, "-"), (0.95, "--")]:
        ax_a.axvline(factors.quantile(q), color="#C44E52", lw=1.0, ls=ls)
    ax_a.set_xlabel("SA(3.0 s) station multiplier")
    ax_a.set_ylabel("Number of stations")
    ax_a.set_title("A Held-out station corrections")
    ax_a.grid(axis="y", color="#e6e6e6")

    colors = {"maximum_negative": "#4C78A8", "q50": "#555555", "maximum_positive": "#C44E52"}
    api_t50 = api_reference[api_reference["period_years"].eq(50)].copy()
    api_t50 = api_t50[np.isfinite(api_t50["annual_rate"].astype(float))]
    for label in ["maximum_negative", "q50", "maximum_positive"]:
        sub = api_t50[api_t50["example_type"].astype(str).eq(label)]
        if sub.empty:
            continue
        ax_b.plot(
            sub["intensity_value"],
            sub["annual_rate"],
            color=colors[label],
            lw=1.4,
            label=f"{label}: {sub['site_code'].iloc[0]}",
        )
    ax_b.set_yscale("log")
    ax_b.set_xlabel("Peak velocity on engineering bedrock (cm/s)")
    ax_b.set_ylabel("Annual exceedance rate")
    ax_b.set_title("B Official J-SHIS hazard-curve API reference")
    ax_b.grid(color="#e6e6e6", which="both")
    ax_b.legend(frameon=False, fontsize=7.4)

    q10 = sa3_values[sa3_values["probability_level"].eq("50y_10pct")].copy()
    ax_c.scatter(
        q10["official_sa3_g"],
        q10["corrected_sa3_g"],
        s=10,
        color="#4C78A8",
        alpha=0.42,
        edgecolors="none",
    )
    limit = float(np.nanmax([q10["official_sa3_g"].max(), q10["corrected_sa3_g"].max()])) * 1.05
    ax_c.plot([0, limit], [0, limit], color="black", lw=0.9)
    ax_c.set_xlim(0, limit)
    ax_c.set_ylim(0, limit)
    ax_c.set_xlabel("Official SA(3.0 s), 50-year 10% (g)")
    ax_c.set_ylabel("Station-corrected SA(3.0 s) (g)")
    ax_c.set_title("C Official response-map ordinates after station correction")
    ax_c.grid(color="#e6e6e6")

    uhs = uhs_examples[
        uhs_examples["probability_level"].eq("50y_10pct")
        & uhs_examples["example_type"].astype(str).isin(["maximum_negative", "q50", "maximum_positive"])
    ].copy()
    line_styles = {"maximum_negative": "-", "q50": "-", "maximum_positive": "-"}
    for label in ["maximum_negative", "q50", "maximum_positive"]:
        sub = uhs[uhs["example_type"].astype(str).eq(label)].sort_values("period_s")
        if sub.empty:
            continue
        ax_d.plot(
            sub["period_s"],
            sub["official_sa_g"],
            color=colors[label],
            lw=1.1,
            ls="--",
            alpha=0.55,
        )
        ax_d.plot(
            sub["period_s"],
            sub["station_corrected_sa_g"],
            color=colors[label],
            lw=1.6,
            ls=line_styles[label],
            label=f"{label}: {sub['site_code'].iloc[0]}",
        )
    ax_d.set_xscale("log")
    ax_d.set_xlim(0.08, 6.0)
    ax_d.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0]))
    ax_d.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax_d.set_xlabel("Period (s)")
    ax_d.set_ylabel("5%-damped SA, 50-year 10% (g)")
    ax_d.set_title("D Official UHS with SA(3.0 s) ordinate adjusted")
    ax_d.grid(color="#e6e6e6", which="both")
    ax_d.legend(frameon=False, fontsize=7.4)

    fig.suptitle("Official-data check of SA(3.0 s) station-correction impacts", fontsize=13.2)
    fig.savefig(figure_dir / "cee_fig8_hazard_impact.pdf", bbox_inches="tight")
    fig.savefig(figure_dir / "cee_fig8_hazard_impact.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--impact", type=Path, default=DEFAULT_IMPACT)
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES)
    parser.add_argument("--site-schema", type=Path, default=DEFAULT_SITE_SCHEMA)
    parser.add_argument("--resp-zip", type=Path, default=DEFAULT_RESP_ZIP)
    parser.add_argument("--api-cache", type=Path, default=DEFAULT_API_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    station_frame = load_station_frame(args.impact, args.site_schema)
    examples = match_examples(args.examples, station_frame)
    api_reference = build_api_reference(examples, args.api_cache)
    meshcodes = set(station_frame["meshcode250_str"].dropna().astype(str))
    response_maps = load_response_maps(args.resp_zip, meshcodes)
    sa3_values = build_sa3_station_values(station_frame, response_maps)
    uhs_examples = build_uhs_examples(examples, response_maps)
    summary = build_summary(sa3_values, examples, api_reference)

    examples.to_csv(args.output_dir / "jshis_official_hazard_response_examples.csv", index=False)
    api_reference.to_csv(args.output_dir / "jshis_official_hazard_curve_api_examples.csv", index=False)
    sa3_values.to_csv(args.output_dir / "jshis_official_response_sa3_station_values.csv", index=False)
    uhs_examples.to_csv(args.output_dir / "jshis_official_response_uhs_examples.csv", index=False)
    summary.to_csv(args.output_dir / "jshis_official_hazard_response_summary.csv", index=False)
    write_report(summary, examples, sa3_values, uhs_examples, args.output_dir / "jshis_official_hazard_response_check.md")
    build_figure(station_frame, examples, api_reference, sa3_values, uhs_examples, args.figure_dir)

    print(f"Wrote {args.output_dir / 'jshis_official_hazard_response_summary.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_official_response_sa3_station_values.csv'}")
    print(f"Wrote {args.output_dir / 'jshis_official_hazard_response_check.md'}")
    print(f"Wrote {args.figure_dir / 'cee_fig8_hazard_impact.pdf'}")


if __name__ == "__main__":
    main()
