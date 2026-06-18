#!/usr/bin/env python3
"""Fetch public NIED K-NET/KiK-net station metadata for current candidates.

The station-list CSV and per-station detail API are public web endpoints.
Bulk waveform HTTPS downloads remain account-based and are intentionally not
handled here.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)
STATIONLIST_PAGE = "https://www.kyoshin.bosai.go.jp/en/stationlist/"
STATIONLIST_DOWNLOAD = (
    "https://www.kyoshin.bosai.go.jp/en/stationlist/download/stationlist/"
)
STATIONINFO_PAGE = "https://www.kyoshin.bosai.go.jp/en/stationinfo/"
STATIONINFO_API = "https://www.kyoshin.bosai.go.jp/en/stationinfo/api/"
USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Safari/537.36"
)


def csrf_from_html(html: str) -> str:
    match = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', html)
    if not match:
        raise RuntimeError("Could not find CSRF token in NIED page.")
    return match.group(1)


def fetch_official_station_list(session: requests.Session) -> pd.DataFrame:
    page = session.get(STATIONLIST_PAGE, timeout=30)
    page.raise_for_status()
    csrf = csrf_from_html(page.text)
    headers = {
        "Referer": STATIONLIST_PAGE,
        "X-CSRFToken": session.cookies.get("csrftoken", csrf),
        "User-Agent": USER_AGENT,
    }
    response = session.post(
        STATIONLIST_DOWNLOAD,
        data={"datakind": "3", "csrfmiddlewaretoken": csrf},
        headers=headers,
        timeout=60,
    )
    response.raise_for_status()
    if "text/csv" not in response.headers.get("content-type", ""):
        raise RuntimeError(
            f"Station list endpoint returned {response.headers.get('content-type')}"
        )
    columns = [
        "network",
        "station_code",
        "station_name",
        "latitude_deg",
        "longitude_deg",
        "elevation_m",
        "depth_m",
        "prefecture",
        "seismograph",
        "notes",
    ]
    frame = pd.read_csv(
        io.StringIO(response.text),
        header=None,
        names=columns,
        quoting=csv.QUOTE_MINIMAL,
    )
    for column in ["latitude_deg", "longitude_deg", "elevation_m", "depth_m"]:
        frame[column] = pd.to_numeric(frame[column].replace("-", pd.NA), errors="coerce")
    return frame


def selected_candidates(output_dir: Path, top_n: int) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    candidate_specs = [
        (
            "knet_high_positive",
            "knet_high_positive_residual_supported.csv",
            "StationCode",
            "K-NET",
            "mean_log10_pga_residual",
            False,
        ),
        (
            "knet_high_negative",
            "knet_high_negative_residual_supported.csv",
            "StationCode",
            "K-NET",
            "mean_log10_pga_residual",
            True,
        ),
        (
            "kiknet_high_ratio",
            "kiknet_high_surface_downhole_ratio.csv",
            "station",
            "KiK-net",
            "horizontal_pga_ratio",
            False,
        ),
        (
            "kiknet_low_ratio",
            "kiknet_low_surface_downhole_ratio.csv",
            "station",
            "KiK-net",
            "horizontal_pga_ratio",
            True,
        ),
    ]
    for group, filename, code_col, network, sort_col, ascending in candidate_specs:
        path = output_dir / filename
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        if sort_col in frame.columns:
            frame = frame.sort_values(sort_col, ascending=ascending)
        frame = frame.head(top_n).copy()
        frame["candidate_group"] = group
        frame["station_code"] = frame[code_col].astype(str)
        frame["network"] = network
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["station_code", "network", "candidate_group"])
    selected = pd.concat(frames, ignore_index=True, sort=False)
    return selected.drop_duplicates(["network", "station_code", "candidate_group"])


def fetch_stationinfo_csrf(session: requests.Session) -> tuple[str, dict[str, str]]:
    page = session.get(STATIONINFO_PAGE, timeout=30)
    page.raise_for_status()
    csrf = csrf_from_html(page.text)
    headers = {
        "Referer": STATIONINFO_PAGE,
        "X-CSRFToken": session.cookies.get("csrftoken", csrf),
        "X-Requested-With": "XMLHttpRequest",
        "User-Agent": USER_AGENT,
    }
    return csrf, headers


def datakind(network: str) -> str:
    if network == "K-NET":
        return "1"
    if network == "KiK-net":
        return "2"
    raise ValueError(f"Unsupported network: {network}")


def fetch_station_detail(
    session: requests.Session,
    csrf: str,
    headers: dict[str, str],
    network: str,
    station_code: str,
) -> dict[str, Any]:
    response = session.post(
        STATIONINFO_API,
        data={
            "datakind": datakind(network),
            "sitecode": station_code,
            "csrfmiddlewaretoken": csrf,
        },
        headers=headers,
        timeout=30,
    )
    response.raise_for_status()
    return response.json()


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip().replace(",", "")
    if value in {"", "-", "-------"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_kik_soil(station_code: str, soil_data: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in soil_data.splitlines():
        match = re.match(
            r"\s*(\d+),\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,\s]+)",
            line,
        )
        if not match:
            continue
        rows.append(
            {
                "station_code": station_code,
                "layer_index": int(match.group(1)),
                "thickness_m": parse_float(match.group(2)),
                "bottom_depth_m": parse_float(match.group(3)),
                "vp_mps": parse_float(match.group(4)),
                "vs_mps": parse_float(match.group(5)),
                "n_value": None,
                "density_g_cm3": None,
                "soil_column": None,
                "profile_format": "kiknet_layer",
            }
        )
    return rows


def parse_knet_soil(station_code: str, soil_data: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    pattern = re.compile(
        r"^\s*(\d+)m\s+"
        r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
        r"(?:\s+(.*?))?\s*$"
    )
    for line in soil_data.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        depth = float(match.group(1))
        rows.append(
            {
                "station_code": station_code,
                "layer_index": int(depth),
                "thickness_m": 1.0,
                "bottom_depth_m": depth,
                "vp_mps": parse_float(match.group(3)),
                "vs_mps": parse_float(match.group(4)),
                "n_value": parse_float(match.group(2)),
                "density_g_cm3": parse_float(match.group(5)),
                "soil_column": (match.group(6) or "").strip() or None,
                "profile_format": "knet_meter",
            }
        )
    return rows


def harmonic_vs_to_depth(profile: pd.DataFrame, depth_m: float) -> float | None:
    if profile.empty:
        return None
    rows = profile.sort_values("bottom_depth_m")
    elapsed = 0.0
    slowness_sum = 0.0
    previous_bottom = 0.0
    last_vs: float | None = None
    for _, row in rows.iterrows():
        vs = row.get("vs_mps")
        bottom = row.get("bottom_depth_m")
        thickness = row.get("thickness_m")
        if pd.isna(vs):
            continue
        if pd.isna(bottom):
            if pd.isna(thickness):
                continue
            bottom = previous_bottom + float(thickness)
        top = previous_bottom
        layer_thickness = max(0.0, min(float(bottom), depth_m) - top)
        if layer_thickness > 0:
            slowness_sum += layer_thickness / float(vs)
            elapsed += layer_thickness
            last_vs = float(vs)
        previous_bottom = max(previous_bottom, float(bottom))
        if elapsed >= depth_m:
            break
    if elapsed < depth_m and last_vs is not None:
        slowness_sum += (depth_m - elapsed) / last_vs
        elapsed = depth_m
    if elapsed <= 0 or slowness_sum <= 0:
        return None
    return depth_m / slowness_sum


def write_report(
    output_dir: Path,
    station_list: pd.DataFrame,
    enriched: pd.DataFrame,
    soil_profiles: pd.DataFrame,
    history: pd.DataFrame,
) -> None:
    knet = station_list[station_list["network"] == "K-NET"]
    kik = station_list[station_list["network"] == "KiK-net"]
    with_soil = enriched[enriched["soil_profile_rows"].fillna(0).gt(0)]
    lines = ["# Official Kyoshin Metadata Acquisition\n\n"]
    lines.append("## What was acquired\n")
    lines.append(
        f"- Official station list: {len(station_list):,} stations "
        f"({len(knet):,} K-NET, {len(kik):,} KiK-net).\n"
    )
    lines.append(
        f"- Candidate station details: {len(enriched):,} station-group rows; "
        f"{with_soil['station_code'].nunique():,} unique stations with parsable soil profiles.\n"
    )
    lines.append(
        f"- Parsed soil-profile rows: {len(soil_profiles):,}; station history rows: {len(history):,}.\n\n"
    )
    lines.append("## Why it matters\n")
    lines.append(
        "The current K-NET/KiK-net candidate stations now have official coordinates, "
        "sensor depths, station names, seismograph types, station histories, and "
        "public per-station soil-profile text where available. This removes the "
        "earlier code-family-only limitation for candidate prioritization.\n\n"
    )
    lines.append("## Candidate highlights\n")
    preview_cols = [
        "candidate_group",
        "network",
        "station_code",
        "station_name",
        "prefecture",
        "latitude_deg",
        "longitude_deg",
        "depth_m",
        "vs30_api_profile_mps",
        "soil_profile_rows",
    ]
    available_cols = [col for col in preview_cols if col in enriched.columns]
    for _, row in enriched[available_cols].head(15).iterrows():
        vs30 = row.get("vs30_api_profile_mps")
        vs30_text = "-" if pd.isna(vs30) else f"{vs30:.0f} m/s"
        depth = row.get("depth_m")
        depth_text = "-" if pd.isna(depth) else f"{depth:.0f} m"
        lines.append(
            f"- {row['candidate_group']} {row['network']} {row['station_code']} "
            f"({row.get('station_name', '-')}, {row.get('prefecture', '-')}): "
            f"lat={row.get('latitude_deg')}, lon={row.get('longitude_deg')}, "
            f"depth={depth_text}, profile VS30 proxy={vs30_text}.\n"
        )
    lines.append("\n## Source and access status\n")
    lines.append(
        "- Station list CSV endpoint: public via the NIED Kyoshin station-list POST endpoint.\n"
    )
    lines.append(
        "- Per-station information API: public with CSRF cookie/session; includes soil-profile text for many stations.\n"
    )
    lines.append(
        "- Bulk waveform HTTPS directories and direct `kyoshin/download/...` paths require NIED user credentials.\n"
    )
    lines.append(
        "- J-SHIS/NIED ground-motion flatfile v2024 is public, but the selected subfile is 255MB and was not fully downloaded in this run because transfer speed was too low.\n\n"
    )
    lines.append("## Research-integrity guardrail\n")
    lines.append(
        "Soil-profile values fetched here are metadata for interpretation and candidate selection. "
        "They do not by themselves prove nonlinear site response or causal amplification. "
        "Claims still need event-wise, leakage-safe validation against waveform-derived or official "
        "ground-motion metrics.\n"
    )
    (output_dir / "official_kyoshin_metadata_acquisition.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--top-n", type=int, default=40)
    parser.add_argument("--sleep-seconds", type=float, default=0.05)
    args = parser.parse_args()

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    station_list = fetch_official_station_list(session)
    station_list.to_csv(args.output_dir / "official_kyoshin_station_list.csv", index=False)

    candidates = selected_candidates(args.output_dir, args.top_n)
    station_lookup = station_list.set_index(["network", "station_code"])

    csrf, headers = fetch_stationinfo_csrf(session)
    enriched_rows: list[dict[str, Any]] = []
    soil_rows: list[dict[str, Any]] = []
    history_rows: list[dict[str, Any]] = []
    seen_detail: set[tuple[str, str]] = set()
    detail_cache: dict[tuple[str, str], dict[str, Any]] = {}

    for _, candidate in candidates.iterrows():
        key = (candidate["network"], candidate["station_code"])
        if key not in seen_detail:
            detail_cache[key] = fetch_station_detail(
                session, csrf, headers, candidate["network"], candidate["station_code"]
            )
            seen_detail.add(key)
            if args.sleep_seconds > 0:
                time.sleep(args.sleep_seconds)
        data = detail_cache[key]
        row = candidate.to_dict()
        if key in station_lookup.index:
            row.update(station_lookup.loc[key].to_dict())
        site_info = data.get("site_info") or {}
        row["station_name_detail"] = site_info.get("sitename")
        row["address"] = site_info.get("address")
        row["soil_image_url"] = (
            f"https://www.kyoshin.bosai.go.jp{data.get('soil_image')}"
            if data.get("soil_image")
            else None
        )
        row["soil_pdf_url"] = (
            f"https://www.kyoshin.bosai.go.jp{data.get('soil_pdf')}"
            if data.get("soil_pdf")
            else None
        )
        soil_data = data.get("soil_data") or ""
        if candidate["network"] == "KiK-net":
            station_soil_rows = parse_kik_soil(candidate["station_code"], soil_data)
        else:
            station_soil_rows = parse_knet_soil(candidate["station_code"], soil_data)
        row["soil_profile_rows"] = len(station_soil_rows)
        if station_soil_rows:
            profile = pd.DataFrame(station_soil_rows)
            row["vs10_api_profile_mps"] = harmonic_vs_to_depth(profile, 10.0)
            row["vs20_api_profile_mps"] = harmonic_vs_to_depth(profile, 20.0)
            row["vs30_api_profile_mps"] = harmonic_vs_to_depth(profile, 30.0)
            for soil_row in station_soil_rows:
                soil_row["network"] = candidate["network"]
                soil_rows.append(soil_row)
        else:
            row["vs10_api_profile_mps"] = None
            row["vs20_api_profile_mps"] = None
            row["vs30_api_profile_mps"] = None

        for history in data.get("site_history", []):
            history_rows.append(
                {
                    "network": candidate["network"],
                    "station_code": candidate["station_code"],
                    "day": history.get("day"),
                    "contents": history.get("contents"),
                }
            )
        enriched_rows.append(row)

    enriched = pd.DataFrame(enriched_rows)
    soil_profiles = pd.DataFrame(soil_rows)
    histories = pd.DataFrame(history_rows)
    enriched.to_csv(
        args.output_dir / "official_candidate_station_metadata.csv", index=False
    )
    soil_profiles.to_csv(
        args.output_dir / "official_candidate_soil_profiles.csv", index=False
    )
    histories.to_csv(
        args.output_dir / "official_candidate_station_history.csv", index=False
    )
    write_report(args.output_dir, station_list, enriched, soil_profiles, histories)


if __name__ == "__main__":
    main()
