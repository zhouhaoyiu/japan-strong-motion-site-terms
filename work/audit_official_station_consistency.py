#!/usr/bin/env python3
"""Audit consistency between local K-NET headers and current NIED metadata."""

from __future__ import annotations

import argparse
import html
import math
import re
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from fetch_official_kyoshin_metadata import (  # noqa: E402
    harmonic_vs_to_depth,
    parse_knet_soil,
)


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)
RELOCATION_PAGE = "https://www.kyoshin.bosai.go.jp/en/stationrelocated/"
SOILINFO_URL = "https://www.kyoshin.bosai.go.jp/en/soilinfo/{station_code}/{variant}/"
USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Safari/537.36"
)


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius_km = 6371.0088
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    )
    return 2.0 * radius_km * math.asin(math.sqrt(a))


def read_relocation_table() -> pd.DataFrame:
    table = pd.read_html(RELOCATION_PAGE)[0]
    table = table.rename(
        columns={
            "Date": "relocation_date",
            "Station Code": "station_code",
            "Station Name": "station_name",
            "Unnamed: 3": "relocation_state",
            "Latitude": "latitude_deg",
            "Longitude": "longitude_deg",
            "Altitude (m)": "elevation_m",
            "Soil condition data (※)": "soil_condition_data_note",
        }
    )
    table["relocation_date"] = pd.to_datetime(
        table["relocation_date"], errors="coerce"
    )
    for column in ["latitude_deg", "longitude_deg", "elevation_m"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["soilinfo_variant"] = table["relocation_state"].map({"Before": 1, "After": 2})
    return table


def summarize_local_knet(output_dir: Path) -> pd.DataFrame:
    header = pd.read_csv(output_dir / "knet_header_index.csv")
    summary = (
        header.groupby("StationCode")
        .agg(
            n_records=("StationCode", "size"),
            min_year=("Year", "min"),
            max_year=("Year", "max"),
            first_origin_time=("OriginTime", "min"),
            last_origin_time=("OriginTime", "max"),
            local_latitude_median=("StationLat", "median"),
            local_longitude_median=("StationLong", "median"),
            local_latitude_min=("StationLat", "min"),
            local_latitude_max=("StationLat", "max"),
            local_longitude_min=("StationLong", "min"),
            local_longitude_max=("StationLong", "max"),
            local_elevation_median=("StationHeight_m", "median"),
            local_elevation_min=("StationHeight_m", "min"),
            local_elevation_max=("StationHeight_m", "max"),
        )
        .reset_index()
        .rename(columns={"StationCode": "station_code"})
    )
    summary["local_coordinate_span_km"] = [
        haversine_km(lat_min, lon_min, lat_max, lon_max)
        for lat_min, lon_min, lat_max, lon_max in zip(
            summary["local_latitude_min"],
            summary["local_longitude_min"],
            summary["local_latitude_max"],
            summary["local_longitude_max"],
        )
    ]
    return summary


def add_relocation_match(
    consistency: pd.DataFrame, relocations: pd.DataFrame, match_threshold_km: float
) -> pd.DataFrame:
    relocation_groups = {
        station: rows.copy() for station, rows in relocations.groupby("station_code")
    }
    rows: list[dict[str, Any]] = []
    for _, row in consistency.iterrows():
        result = row.to_dict()
        station_relocations = relocation_groups.get(row["station_code"])
        result["has_official_relocation_entry"] = station_relocations is not None
        result["relocation_dates"] = None
        result["closest_relocation_state"] = None
        result["closest_relocation_date"] = None
        result["closest_relocation_distance_km"] = None
        result["closest_relocation_latitude_deg"] = None
        result["closest_relocation_longitude_deg"] = None
        result["closest_soilinfo_variant"] = None
        if station_relocations is not None and not station_relocations.empty:
            result["relocation_dates"] = ";".join(
                d.strftime("%Y-%m-%d")
                for d in station_relocations["relocation_date"].dropna().drop_duplicates()
            )
            distances = []
            for _, rel in station_relocations.iterrows():
                distance = haversine_km(
                    row["local_latitude_median"],
                    row["local_longitude_median"],
                    rel["latitude_deg"],
                    rel["longitude_deg"],
                )
                distances.append((distance, rel))
            distance, rel = min(distances, key=lambda item: item[0])
            result["closest_relocation_state"] = rel["relocation_state"]
            result["closest_relocation_date"] = rel["relocation_date"]
            result["closest_relocation_distance_km"] = distance
            result["closest_relocation_latitude_deg"] = rel["latitude_deg"]
            result["closest_relocation_longitude_deg"] = rel["longitude_deg"]
            result["closest_soilinfo_variant"] = rel["soilinfo_variant"]
        current_delta = result["current_coordinate_delta_km"]
        closest_delta = result["closest_relocation_distance_km"]
        if pd.isna(current_delta):
            recommendation = "missing_current_official_metadata"
        elif current_delta <= match_threshold_km:
            recommendation = "current_metadata_matches_local_header"
        elif closest_delta is not None and closest_delta <= match_threshold_km:
            recommendation = (
                f"use_{str(result['closest_relocation_state']).lower()}_relocation_metadata"
            )
        elif current_delta >= 1.0:
            recommendation = "high_mismatch_do_not_use_current_site_metadata_directly"
        else:
            recommendation = "moderate_mismatch_verify_before_site_parameter_join"
        result["metadata_join_recommendation"] = recommendation
        rows.append(result)
    return pd.DataFrame(rows)


def build_consistency(output_dir: Path, match_threshold_km: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    local = summarize_local_knet(output_dir)
    official = pd.read_csv(output_dir / "official_kyoshin_station_list.csv")
    official = official[official["network"] == "K-NET"].copy()
    official = official.rename(
        columns={
            "latitude_deg": "current_official_latitude_deg",
            "longitude_deg": "current_official_longitude_deg",
            "elevation_m": "current_official_elevation_m",
            "station_name": "current_official_station_name",
            "prefecture": "current_official_prefecture",
            "seismograph": "current_official_seismograph",
        }
    )
    consistency = local.merge(
        official[
            [
                "station_code",
                "current_official_station_name",
                "current_official_latitude_deg",
                "current_official_longitude_deg",
                "current_official_elevation_m",
                "current_official_prefecture",
                "current_official_seismograph",
            ]
        ],
        on="station_code",
        how="left",
    )
    consistency["current_coordinate_delta_km"] = [
        (
            haversine_km(lat, lon, current_lat, current_lon)
            if not pd.isna(current_lat)
            else pd.NA
        )
        for lat, lon, current_lat, current_lon in zip(
            consistency["local_latitude_median"],
            consistency["local_longitude_median"],
            consistency["current_official_latitude_deg"],
            consistency["current_official_longitude_deg"],
        )
    ]
    consistency["current_elevation_delta_m"] = (
        consistency["current_official_elevation_m"]
        - consistency["local_elevation_median"]
    )
    relocations = read_relocation_table()
    consistency = add_relocation_match(consistency, relocations, match_threshold_km)
    return consistency, relocations


def extract_tag_text(html_text: str, tag_id: str) -> str | None:
    pattern = re.compile(rf'id="{re.escape(tag_id)}"[^>]*>(.*?)</', re.DOTALL)
    match = pattern.search(html_text)
    if not match:
        return None
    return html.unescape(re.sub(r"<[^>]+>", "", match.group(1))).strip()


def fetch_historical_soilinfo(
    session: requests.Session,
    station_code: str,
    variant: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    url = SOILINFO_URL.format(station_code=station_code, variant=variant)
    response = session.get(url, timeout=30)
    response.raise_for_status()
    text = response.text
    pre_match = re.search(
        r'<pre id="soil_quality_data">(.*?)</pre>', text, flags=re.DOTALL
    )
    soil_data = html.unescape(pre_match.group(1)) if pre_match else ""
    soil_data = re.sub(r"<[^>]+>", "", soil_data)
    rows = parse_knet_soil(station_code, soil_data)
    detail = {
        "station_code": station_code,
        "soilinfo_variant": variant,
        "soilinfo_url": url,
        "soilinfo_station_name": extract_tag_text(text, "sitename"),
        "soilinfo_latitude_deg": extract_tag_text(text, "lat"),
        "soilinfo_longitude_deg": extract_tag_text(text, "lon"),
        "soilinfo_elevation_m": extract_tag_text(text, "elevation"),
        "soilinfo_num_relocation": extract_tag_text(text, "num_relocation"),
        "soil_profile_rows": len(rows),
    }
    profile = pd.DataFrame(rows)
    detail["vs10_soilinfo_profile_mps"] = harmonic_vs_to_depth(profile, 10.0)
    detail["vs20_soilinfo_profile_mps"] = harmonic_vs_to_depth(profile, 20.0)
    detail["vs30_soilinfo_profile_mps"] = harmonic_vs_to_depth(profile, 30.0)
    for row in rows:
        row["soilinfo_variant"] = variant
        row["soilinfo_url"] = url
    return detail, rows


def fetch_candidate_historical_soils(
    output_dir: Path,
    consistency: pd.DataFrame,
    match_threshold_km: float,
    sleep_seconds: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(output_dir / "official_candidate_station_metadata.csv")
    knet_candidates = candidates[candidates["network"] == "K-NET"][
        ["station_code", "candidate_group"]
    ].drop_duplicates()
    candidate_consistency = knet_candidates.merge(
        consistency,
        on="station_code",
        how="left",
    )
    to_fetch = candidate_consistency[
        candidate_consistency["closest_relocation_distance_km"].notna()
        & candidate_consistency["closest_relocation_distance_km"].le(match_threshold_km)
        & candidate_consistency["closest_soilinfo_variant"].notna()
        & (
            candidate_consistency["metadata_join_recommendation"].str.startswith(
                "use_"
            )
        )
    ].copy()
    unique_fetch = to_fetch[
        ["station_code", "closest_soilinfo_variant", "closest_relocation_state"]
    ].drop_duplicates()
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    summaries: list[dict[str, Any]] = []
    profiles: list[dict[str, Any]] = []
    for _, row in unique_fetch.iterrows():
        detail, rows = fetch_historical_soilinfo(
            session, row["station_code"], int(row["closest_soilinfo_variant"])
        )
        detail["matched_relocation_state"] = row["closest_relocation_state"]
        summaries.append(detail)
        profiles.extend(rows)
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)
    return pd.DataFrame(summaries), pd.DataFrame(profiles)


def write_report(
    output_dir: Path,
    consistency: pd.DataFrame,
    relocations: pd.DataFrame,
    candidate_consistency: pd.DataFrame,
    historical_soil_summary: pd.DataFrame,
    match_threshold_km: float,
) -> None:
    current_delta = consistency["current_coordinate_delta_km"]
    over_threshold = current_delta.gt(match_threshold_km)
    over_1km = current_delta.gt(1.0)
    with_relocation = consistency["has_official_relocation_entry"].fillna(False)
    lines = ["# Station Metadata Consistency Audit\n\n"]
    lines.append("## Main finding\n")
    lines.append(
        f"The local K-NET records span {int(consistency['min_year'].min())}-"
        f"{int(consistency['max_year'].max())}, while the official station-list "
        "endpoint returns current station metadata. These are not always the same "
        "physical site.\n\n"
    )
    lines.append("## Coordinate mismatch statistics\n")
    lines.append(f"- Local K-NET stations audited: {len(consistency):,}.\n")
    lines.append(f"- Official relocation entries parsed: {len(relocations) // 2:,} relocations.\n")
    lines.append(
        f"- Current official coordinate mismatch > {match_threshold_km:.1f} km: "
        f"{int(over_threshold.sum()):,} stations.\n"
    )
    lines.append(f"- Current official coordinate mismatch > 1.0 km: {int(over_1km.sum()):,} stations.\n")
    lines.append(
        f"- Mismatch > {match_threshold_km:.1f} km with official relocation entry: "
        f"{int((over_threshold & with_relocation).sum()):,} stations.\n"
    )
    lines.append(
        f"- Mismatch > {match_threshold_km:.1f} km without relocation entry on this page: "
        f"{int((over_threshold & ~with_relocation).sum()):,} stations.\n\n"
    )
    lines.append("## Largest local-vs-current coordinate differences\n")
    largest = consistency.sort_values("current_coordinate_delta_km", ascending=False).head(12)
    for _, row in largest.iterrows():
        lines.append(
            f"- {row['station_code']} ({row['current_official_station_name']}): "
            f"delta={row['current_coordinate_delta_km']:.2f} km, "
            f"local=({row['local_latitude_median']:.4f}, {row['local_longitude_median']:.4f}), "
            f"current=({row['current_official_latitude_deg']:.4f}, "
            f"{row['current_official_longitude_deg']:.4f}), "
            f"recommendation={row['metadata_join_recommendation']}.\n"
        )
    lines.append("\n## Candidate station implications\n")
    candidate_bad = candidate_consistency[
        candidate_consistency["current_coordinate_delta_km"].gt(match_threshold_km)
    ].sort_values("current_coordinate_delta_km", ascending=False)
    lines.append(
        f"- K-NET candidate station rows checked: {len(candidate_consistency):,}; "
        f"{len(candidate_bad):,} exceed {match_threshold_km:.1f} km against current metadata.\n"
    )
    if not historical_soil_summary.empty:
        lines.append(
            f"- Historical relocation-state soil profiles fetched for "
            f"{historical_soil_summary['station_code'].nunique():,} candidate stations.\n"
        )
    for _, row in candidate_bad.head(12).iterrows():
        lines.append(
            f"- {row['candidate_group']} {row['station_code']}: "
            f"current delta={row['current_coordinate_delta_km']:.2f} km, "
            f"closest relocation state={row['closest_relocation_state']} "
            f"({row['closest_relocation_distance_km']:.2f} km), "
            f"recommendation={row['metadata_join_recommendation']}.\n"
        )
    lines.append("\n## Methodological decision\n")
    lines.append(
        "For the 1997-2006 local K-NET waveform dataset, model geometry should continue "
        "to use per-record header coordinates. Current official station metadata and "
        "soil profiles should be joined only when the current coordinates match the "
        "header coordinates, or when a matched before/after relocation-state profile "
        "is explicitly used. This is now enforced as a metadata guardrail rather than "
        "a cosmetic data-cleaning issue.\n\n"
    )
    lines.append("## Research-integrity guardrail\n")
    lines.append(
        "A station code is not sufficient evidence that two records share the same physical "
        "site across time. Relocations, coordinate changes, and station-history rows must "
        "be handled before interpreting station terms as soil/site effects.\n"
    )
    (output_dir / "station_metadata_consistency_audit.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--match-threshold-km", type=float, default=0.2)
    parser.add_argument("--sleep-seconds", type=float, default=0.05)
    args = parser.parse_args()

    consistency, relocations = build_consistency(args.output_dir, args.match_threshold_km)
    relocations.to_csv(args.output_dir / "official_kyoshin_relocated_station_list.csv", index=False)
    consistency.to_csv(args.output_dir / "official_knet_station_coordinate_consistency.csv", index=False)

    candidates = pd.read_csv(args.output_dir / "official_candidate_station_metadata.csv")
    knet_candidates = candidates[candidates["network"] == "K-NET"][
        ["station_code", "candidate_group"]
    ].drop_duplicates()
    candidate_consistency = knet_candidates.merge(
        consistency,
        on="station_code",
        how="left",
    )
    candidate_consistency.to_csv(
        args.output_dir / "official_candidate_knet_metadata_consistency.csv",
        index=False,
    )

    historical_soil_summary, historical_soil_profiles = fetch_candidate_historical_soils(
        args.output_dir,
        consistency,
        args.match_threshold_km,
        args.sleep_seconds,
    )
    historical_soil_summary.to_csv(
        args.output_dir / "official_candidate_relocation_soilinfo_summary.csv",
        index=False,
    )
    historical_soil_profiles.to_csv(
        args.output_dir / "official_candidate_relocation_soil_profiles.csv",
        index=False,
    )

    write_report(
        args.output_dir,
        consistency,
        relocations,
        candidate_consistency,
        historical_soil_summary,
        args.match_threshold_km,
    )


if __name__ == "__main__":
    main()
