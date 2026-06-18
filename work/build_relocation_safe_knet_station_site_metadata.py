#!/usr/bin/env python3
"""Build relocation-safe K-NET site metadata for station-term interpretation."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

import sys

sys.path.insert(0, str(Path(__file__).parent))
from audit_official_station_consistency import fetch_historical_soilinfo  # noqa: E402
from fetch_official_kyoshin_metadata import (  # noqa: E402
    STATIONINFO_API,
    datakind,
    fetch_stationinfo_csrf,
    harmonic_vs_to_depth,
    parse_knet_soil,
)


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)
DEFAULT_CACHE_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/work/external_data/kyoshin_site_cache"
)
SAFE_RECOMMENDATIONS = {
    "current_metadata_matches_local_header",
    "use_before_relocation_metadata",
    "use_after_relocation_metadata",
}


def fetch_current_station_detail(
    session: requests.Session,
    csrf: str,
    headers: dict[str, str],
    cache_dir: Path,
    station_code: str,
    use_cache: bool,
) -> dict[str, Any]:
    cache_path = cache_dir / f"current_knet_{station_code}.json"
    if use_cache and cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    response = session.post(
        STATIONINFO_API,
        data={
            "datakind": datakind("K-NET"),
            "sitecode": station_code,
            "csrfmiddlewaretoken": csrf,
        },
        headers=headers,
        timeout=30,
    )
    response.raise_for_status()
    data = response.json()
    cache_path.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    return data


def current_detail_to_profile(
    station_code: str,
    detail: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = parse_knet_soil(station_code, detail.get("soil_data") or "")
    profile = pd.DataFrame(rows)
    site_info = detail.get("site_info") or {}
    summary = {
        "station_code": station_code,
        "soil_profile_source": "current_stationinfo_coordinate_matched",
        "soil_profile_reference_url": (
            f"https://www.kyoshin.bosai.go.jp{detail.get('soil_pdf')}"
            if detail.get("soil_pdf")
            else None
        ),
        "soil_profile_rows": len(rows),
        "safe_site_latitude_deg": pd.to_numeric(site_info.get("lat"), errors="coerce"),
        "safe_site_longitude_deg": pd.to_numeric(site_info.get("lon"), errors="coerce"),
        "safe_site_elevation_m": pd.to_numeric(site_info.get("elevation"), errors="coerce"),
        "safe_vs10_proxy_mps": harmonic_vs_to_depth(profile, 10.0),
        "safe_vs20_proxy_mps": harmonic_vs_to_depth(profile, 20.0),
        "safe_vs30_proxy_mps": harmonic_vs_to_depth(profile, 30.0),
    }
    for row in rows:
        row["safe_profile_source"] = "current_stationinfo_coordinate_matched"
        row["soil_profile_reference_url"] = summary["soil_profile_reference_url"]
    return summary, rows


def relocation_detail_to_profile(
    station_code: str,
    state: str,
    variant: int,
    detail: dict[str, Any],
    rows: list[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source = f"relocation_{state.lower()}_soilinfo_coordinate_matched"
    summary = {
        "station_code": station_code,
        "soil_profile_source": source,
        "soil_profile_reference_url": detail["soilinfo_url"],
        "soil_profile_rows": detail["soil_profile_rows"],
        "safe_site_latitude_deg": pd.to_numeric(
            detail["soilinfo_latitude_deg"], errors="coerce"
        ),
        "safe_site_longitude_deg": pd.to_numeric(
            detail["soilinfo_longitude_deg"], errors="coerce"
        ),
        "safe_site_elevation_m": pd.to_numeric(
            detail["soilinfo_elevation_m"], errors="coerce"
        ),
        "safe_vs10_proxy_mps": detail["vs10_soilinfo_profile_mps"],
        "safe_vs20_proxy_mps": detail["vs20_soilinfo_profile_mps"],
        "safe_vs30_proxy_mps": detail["vs30_soilinfo_profile_mps"],
        "soilinfo_variant": variant,
    }
    for row in rows:
        row["safe_profile_source"] = source
        row["soil_profile_reference_url"] = detail["soilinfo_url"]
    return summary, rows


def build_target_table(output_dir: Path, min_records: int) -> pd.DataFrame:
    consistency = pd.read_csv(output_dir / "official_knet_station_coordinate_consistency.csv")
    adjusted = pd.read_csv(output_dir / "knet_adjusted_station_terms.csv").rename(
        columns={"n_records": "adjusted_term_n_records"}
    )
    residual = pd.read_csv(output_dir / "knet_station_residual_terms.csv").rename(
        columns={"n_records": "residual_n_records"}
    )
    targets = consistency.merge(
        adjusted,
        left_on="station_code",
        right_on="StationCode",
        how="inner",
    ).merge(
        residual,
        left_on="station_code",
        right_on="StationCode",
        how="left",
        suffixes=("", "_residual"),
    )
    targets = targets[
        targets["adjusted_term_n_records"].ge(min_records)
        & targets["metadata_join_recommendation"].isin(SAFE_RECOMMENDATIONS)
    ].copy()
    return targets.sort_values(["metadata_join_recommendation", "station_code"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--min-records", type=int, default=10)
    parser.add_argument("--sleep-seconds", type=float, default=0.05)
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args()

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    targets = build_target_table(args.output_dir, args.min_records)

    session = requests.Session()
    csrf, headers = fetch_stationinfo_csrf(session)
    summaries: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []

    use_cache = not args.no_cache
    for idx, row in targets.iterrows():
        station_code = row["station_code"]
        recommendation = row["metadata_join_recommendation"]
        if recommendation == "current_metadata_matches_local_header":
            detail = fetch_current_station_detail(
                session, csrf, headers, args.cache_dir, station_code, use_cache
            )
            summary, rows = current_detail_to_profile(station_code, detail)
        else:
            variant = int(row["closest_soilinfo_variant"])
            state = str(row["closest_relocation_state"])
            detail, rows = fetch_historical_soilinfo(session, station_code, variant)
            summary, rows = relocation_detail_to_profile(
                station_code, state, variant, detail, rows
            )
        summary.update(
            {
                "metadata_join_recommendation": recommendation,
                "regularized_term_log10_pga": row["regularized_term_log10_pga"],
                "mean_log10_pga_residual": row["mean_log10_pga_residual"],
                "median_log10_pga_residual": row["median_log10_pga_residual"],
                "adjusted_term_n_records": row["adjusted_term_n_records"],
                "median_distance_km": row["median_distance_km"],
                "median_magnitude": row["median_magnitude"],
                "current_coordinate_delta_km": row["current_coordinate_delta_km"],
                "closest_relocation_distance_km": row["closest_relocation_distance_km"],
                "local_latitude_median": row["local_latitude_median"],
                "local_longitude_median": row["local_longitude_median"],
                "local_elevation_median": row["local_elevation_median"],
            }
        )
        summary["safe_for_site_parameter_join"] = (
            summary["soil_profile_rows"] > 0
            and not pd.isna(summary["safe_vs30_proxy_mps"])
        )
        summaries.append(summary)
        for profile_row in rows:
            profile_row["station_code"] = station_code
            profile_rows.append(profile_row)

        if (idx + 1) % 25 == 0 or idx == len(targets) - 1:
            print(f"Processed {idx + 1}/{len(targets)} station profiles", flush=True)
        if args.sleep_seconds > 0:
            time.sleep(args.sleep_seconds)

    summary_df = pd.DataFrame(summaries)
    profile_df = pd.DataFrame(profile_rows)
    summary_path = (
        args.output_dir
        / f"relocation_safe_knet_station_site_parameters_min{args.min_records}.csv"
    )
    profile_path = (
        args.output_dir
        / f"relocation_safe_knet_station_soil_profiles_min{args.min_records}.csv"
    )
    summary_df.to_csv(summary_path, index=False)
    profile_df.to_csv(profile_path, index=False)

    safe_count = int(summary_df["safe_for_site_parameter_join"].fillna(False).sum())
    report = [
        "# Relocation-Safe K-NET Station Site Metadata\n\n",
        f"- Minimum station-term record count: {args.min_records}.\n",
        f"- Coordinate/relocation-safe stations requested: {len(targets):,}.\n",
        f"- Stations with parsable soil profile and VS30 proxy: {safe_count:,}.\n",
        f"- Soil-profile rows parsed: {len(profile_df):,}.\n",
        "\n## Source counts\n",
    ]
    for source, count in summary_df["soil_profile_source"].value_counts().items():
        report.append(f"- {source}: {int(count):,}.\n")
    missing = summary_df[~summary_df["safe_for_site_parameter_join"].fillna(False)]
    if not missing.empty:
        report.append("\n## Safe coordinate but missing soil profile\n")
        for _, row in missing.iterrows():
            report.append(f"- {row['station_code']}: {row['soil_profile_source']}.\n")
    report.append("\n## Guardrail\n")
    report.append(
        "These profile-derived VS proxies are used for station-term interpretation only. "
        "Final manuscript-scale modeling should prefer official J-SHIS/NIED flatfile "
        "`site_schema.tsv` values when available.\n"
    )
    (
        args.output_dir
        / f"relocation_safe_knet_station_site_metadata_min{args.min_records}.md"
    ).write_text("".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
