#!/usr/bin/env python3
"""Build relocation-safe site metadata tables for candidate stations."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)


def first_existing(frame: pd.DataFrame, columns: list[str]) -> str | None:
    for column in columns:
        if column in frame.columns:
            return column
    return None


def assign_knet_safe_metadata(output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(output_dir / "official_candidate_station_metadata.csv")
    consistency = pd.read_csv(output_dir / "official_candidate_knet_metadata_consistency.csv")
    relocation_soil = pd.read_csv(output_dir / "official_candidate_relocation_soilinfo_summary.csv")
    current_profiles = pd.read_csv(output_dir / "official_candidate_soil_profiles.csv")
    relocation_profiles = pd.read_csv(output_dir / "official_candidate_relocation_soil_profiles.csv")

    knet = candidates[candidates["network"] == "K-NET"].copy()
    kik = candidates[candidates["network"] == "KiK-net"].copy()

    consistency_cols = [
        "station_code",
        "current_coordinate_delta_km",
        "current_elevation_delta_m",
        "has_official_relocation_entry",
        "relocation_dates",
        "closest_relocation_state",
        "closest_relocation_date",
        "closest_relocation_distance_km",
        "closest_soilinfo_variant",
        "metadata_join_recommendation",
    ]
    knet = knet.merge(
        consistency[consistency_cols],
        on="station_code",
        how="left",
    )
    relocation_cols = [
        "station_code",
        "soilinfo_variant",
        "soilinfo_url",
        "soilinfo_latitude_deg",
        "soilinfo_longitude_deg",
        "soilinfo_elevation_m",
        "soil_profile_rows",
        "vs10_soilinfo_profile_mps",
        "vs20_soilinfo_profile_mps",
        "vs30_soilinfo_profile_mps",
        "matched_relocation_state",
    ]
    relocation_soil = relocation_soil[relocation_cols].rename(
        columns={
            "soil_profile_rows": "relocation_soil_profile_rows",
            "soilinfo_url": "relocation_soilinfo_url",
        }
    )
    knet = knet.merge(
        relocation_soil,
        left_on=["station_code", "closest_soilinfo_variant"],
        right_on=["station_code", "soilinfo_variant"],
        how="left",
    )

    safe_rows = []
    for _, row in knet.iterrows():
        recommendation = row.get("metadata_join_recommendation")
        safe = row.to_dict()
        safe["safe_for_coordinate_join"] = False
        safe["safe_for_site_parameter_join"] = False
        safe["site_parameter_source"] = "unsafe_not_joined"
        safe["safe_site_latitude_deg"] = row.get("station_latitude")
        safe["safe_site_longitude_deg"] = row.get("station_longitude")
        safe["safe_site_elevation_m"] = row.get("station_height_m")
        safe["safe_vs10_proxy_mps"] = pd.NA
        safe["safe_vs20_proxy_mps"] = pd.NA
        safe["safe_vs30_proxy_mps"] = pd.NA
        safe["safe_soil_profile_rows"] = 0
        safe["safe_soil_profile_reference_url"] = pd.NA

        if recommendation == "current_metadata_matches_local_header":
            safe["safe_for_coordinate_join"] = True
            safe["safe_site_latitude_deg"] = row.get("latitude_deg")
            safe["safe_site_longitude_deg"] = row.get("longitude_deg")
            safe["safe_site_elevation_m"] = row.get("elevation_m")
            has_profile = (
                pd.to_numeric(row.get("soil_profile_rows"), errors="coerce") > 0
                and not pd.isna(row.get("vs30_api_profile_mps"))
            )
            if has_profile:
                safe["safe_for_site_parameter_join"] = True
                safe["site_parameter_source"] = "current_stationinfo_coordinate_matched"
                safe["safe_vs10_proxy_mps"] = row.get("vs10_api_profile_mps")
                safe["safe_vs20_proxy_mps"] = row.get("vs20_api_profile_mps")
                safe["safe_vs30_proxy_mps"] = row.get("vs30_api_profile_mps")
                safe["safe_soil_profile_rows"] = row.get("soil_profile_rows")
                safe["safe_soil_profile_reference_url"] = row.get("soil_pdf_url")
            else:
                safe["site_parameter_source"] = (
                    "current_coordinate_matched_but_soil_profile_missing"
                )
        elif isinstance(recommendation, str) and recommendation.startswith("use_"):
            has_relocation_profile = not pd.isna(row.get("vs30_soilinfo_profile_mps"))
            if has_relocation_profile:
                state = str(row.get("matched_relocation_state")).lower()
                safe["safe_for_coordinate_join"] = True
                safe["safe_for_site_parameter_join"] = True
                safe["site_parameter_source"] = f"relocation_{state}_soilinfo_coordinate_matched"
                safe["safe_site_latitude_deg"] = pd.to_numeric(
                    row.get("soilinfo_latitude_deg"), errors="coerce"
                )
                safe["safe_site_longitude_deg"] = pd.to_numeric(
                    row.get("soilinfo_longitude_deg"), errors="coerce"
                )
                safe["safe_site_elevation_m"] = pd.to_numeric(
                    row.get("soilinfo_elevation_m"), errors="coerce"
                )
                safe["safe_vs10_proxy_mps"] = row.get("vs10_soilinfo_profile_mps")
                safe["safe_vs20_proxy_mps"] = row.get("vs20_soilinfo_profile_mps")
                safe["safe_vs30_proxy_mps"] = row.get("vs30_soilinfo_profile_mps")
                safe["safe_soil_profile_rows"] = row.get("relocation_soil_profile_rows")
                safe["safe_soil_profile_reference_url"] = row.get(
                    "relocation_soilinfo_url"
                )
            else:
                safe["site_parameter_source"] = "relocation_match_but_profile_missing"
        elif isinstance(recommendation, str) and recommendation.startswith("moderate"):
            safe["site_parameter_source"] = "unsafe_moderate_coordinate_mismatch_verify"
        elif isinstance(recommendation, str) and recommendation.startswith("high"):
            safe["site_parameter_source"] = "unsafe_high_coordinate_mismatch"

        safe_rows.append(safe)

    safe_knet = pd.DataFrame(safe_rows)
    safe_knet.to_csv(output_dir / "relocation_safe_knet_candidate_site_parameters.csv", index=False)

    # Build profile table with only the profiles allowed by the safe join logic.
    safe_current_codes = set(
        safe_knet.loc[
            safe_knet["site_parameter_source"].eq("current_stationinfo_coordinate_matched"),
            "station_code",
        ]
    )
    current_safe_profiles = current_profiles[
        (current_profiles["network"] == "K-NET")
        & current_profiles["station_code"].isin(safe_current_codes)
    ].copy()
    current_safe_profiles["safe_profile_source"] = "current_stationinfo_coordinate_matched"

    safe_relocation = safe_knet[
        safe_knet["site_parameter_source"].str.startswith("relocation_", na=False)
    ][["station_code", "closest_soilinfo_variant", "site_parameter_source"]].drop_duplicates()
    relocation_safe_profiles = relocation_profiles.merge(
        safe_relocation,
        left_on=["station_code", "soilinfo_variant"],
        right_on=["station_code", "closest_soilinfo_variant"],
        how="inner",
    )
    relocation_safe_profiles = relocation_safe_profiles.drop(
        columns=["closest_soilinfo_variant"], errors="ignore"
    ).rename(columns={"site_parameter_source": "safe_profile_source"})
    safe_profiles = pd.concat(
        [current_safe_profiles, relocation_safe_profiles],
        ignore_index=True,
        sort=False,
    )
    safe_profiles.to_csv(output_dir / "relocation_safe_knet_candidate_soil_profiles.csv", index=False)

    # KiK-net candidates lack local station coordinates in the converted one-event miniSEED set.
    # Keep them in a combined table but mark the status explicitly.
    if not kik.empty:
        kik = kik.copy()
        kik["safe_for_site_parameter_join"] = pd.NA
        kik["site_parameter_source"] = "kiknet_current_metadata_not_relocation_audited"
        kik["safe_site_latitude_deg"] = kik["latitude_deg"]
        kik["safe_site_longitude_deg"] = kik["longitude_deg"]
        kik["safe_site_elevation_m"] = kik["elevation_m"]
        kik["safe_vs10_proxy_mps"] = kik["vs10_api_profile_mps"]
        kik["safe_vs20_proxy_mps"] = kik["vs20_api_profile_mps"]
        kik["safe_vs30_proxy_mps"] = kik["vs30_api_profile_mps"]
        kik["safe_soil_profile_rows"] = kik["soil_profile_rows"]
        kik["safe_soil_profile_reference_url"] = kik["soil_pdf_url"]
    combined = pd.concat([safe_knet, kik], ignore_index=True, sort=False)
    combined.to_csv(output_dir / "relocation_safe_candidate_site_metadata.csv", index=False)
    return safe_knet, safe_profiles


def write_report(output_dir: Path, safe_knet: pd.DataFrame, safe_profiles: pd.DataFrame) -> None:
    source_counts = safe_knet["site_parameter_source"].value_counts(dropna=False)
    safe_count = int(safe_knet["safe_for_site_parameter_join"].fillna(False).sum())
    unsafe = safe_knet[~safe_knet["safe_for_site_parameter_join"].fillna(False)].copy()
    relocation_fixed = safe_knet[
        safe_knet["site_parameter_source"].str.startswith("relocation_", na=False)
    ].copy()
    lines = ["# Relocation-Safe Site Metadata\n\n"]
    lines.append("## Decision implemented\n")
    lines.append(
        "For the 1997-2006 K-NET waveform records, station-code joins are no longer "
        "accepted by default. Site parameters are allowed only when official current "
        "coordinates match the local header coordinates, or when a matching official "
        "before/after relocation soil profile is used.\n\n"
    )
    lines.append("## K-NET candidate status\n")
    lines.append(f"- K-NET candidate rows: {len(safe_knet):,}.\n")
    lines.append(f"- Safe for site-parameter join: {safe_count:,}.\n")
    lines.append(f"- Unsafe pending manual verification: {len(unsafe):,}.\n")
    lines.append(f"- Safe soil-profile rows retained: {len(safe_profiles):,}.\n")
    for source, count in source_counts.items():
        lines.append(f"- {source}: {int(count):,}.\n")
    lines.append("\n## Relocation corrections with large VS30 changes\n")
    if not relocation_fixed.empty:
        relocation_fixed["current_minus_safe_vs30_proxy_mps"] = (
            pd.to_numeric(relocation_fixed["vs30_api_profile_mps"], errors="coerce")
            - pd.to_numeric(relocation_fixed["safe_vs30_proxy_mps"], errors="coerce")
        )
        preview = relocation_fixed.reindex(
            relocation_fixed["current_minus_safe_vs30_proxy_mps"].abs().sort_values(
                ascending=False
            ).index
        ).head(10)
        for _, row in preview.iterrows():
            lines.append(
                f"- {row['station_code']}: current VS30 proxy="
                f"{row['vs30_api_profile_mps']:.0f} m/s, "
                f"relocation-safe VS30 proxy={row['safe_vs30_proxy_mps']:.0f} m/s, "
                f"delta={row['current_minus_safe_vs30_proxy_mps']:.0f} m/s, "
                f"source={row['site_parameter_source']}.\n"
            )
    lines.append("\n## Unsafe K-NET candidates not joined\n")
    for _, row in unsafe.sort_values("current_coordinate_delta_km", ascending=False).iterrows():
        lines.append(
            f"- {row['station_code']}: current coordinate delta="
            f"{row['current_coordinate_delta_km']:.2f} km; "
            f"source={row['site_parameter_source']}.\n"
        )
    lines.append("\n## Output files\n")
    lines.append("- `relocation_safe_knet_candidate_site_parameters.csv`\n")
    lines.append("- `relocation_safe_knet_candidate_soil_profiles.csv`\n")
    lines.append("- `relocation_safe_candidate_site_metadata.csv`\n\n")
    lines.append("## Figures\n")
    lines.append("- `figures/knet_station_coordinate_mismatch_histogram.png`\n")
    lines.append("- `figures/knet_candidate_site_metadata_source_counts.png`\n")
    lines.append("- `figures/knet_relocation_vs30_correction.png`\n\n")
    lines.append("## Guardrail\n")
    lines.append(
        "The VS10/VS20/VS30 values in these tables are profile-derived proxies from "
        "the public Kyoshin station/soil pages. They are suitable for candidate "
        "screening and interpretation, not a substitute for the official J-SHIS "
        "flatfile site_schema.tsv values once that dataset is downloaded.\n"
    )
    (output_dir / "relocation_safe_site_metadata.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def write_figures(output_dir: Path, safe_knet: pd.DataFrame) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)

    consistency_path = output_dir / "official_knet_station_coordinate_consistency.csv"
    if consistency_path.exists():
        consistency = pd.read_csv(consistency_path)
        delta = pd.to_numeric(
            consistency["current_coordinate_delta_km"], errors="coerce"
        ).dropna()
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        ax.hist(
            delta,
            bins=[0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 15],
            color="#4267ac",
        )
        ax.axvline(
            0.2,
            color="#b13f3f",
            linestyle="--",
            linewidth=1.2,
            label="0.2 km threshold",
        )
        ax.set_xscale("log")
        ax.set_xlabel("Distance between local header and current coordinates (km)")
        ax.set_ylabel("Number of K-NET stations")
        ax.set_title("K-NET Station Coordinate Consistency")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(
            figure_dir / "knet_station_coordinate_mismatch_histogram.png",
            dpi=180,
        )
        plt.close(fig)

    counts = safe_knet["site_parameter_source"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(7.8, 3.9))
    colors = [
        "#7c7c7c" if ("unsafe" in str(idx) or "missing" in str(idx)) else "#2f7d5c"
        for idx in counts.index
    ]
    ax.barh(counts.index, counts.values, color=colors)
    ax.set_xlabel("K-NET candidate rows")
    ax.set_title("Relocation-Safe Metadata Source")
    for y, value in enumerate(counts.values):
        ax.text(value + 0.3, y, str(int(value)), va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(
        figure_dir / "knet_candidate_site_metadata_source_counts.png",
        dpi=180,
    )
    plt.close(fig)

    relocation_fixed = safe_knet[
        safe_knet["site_parameter_source"].str.startswith("relocation_", na=False)
    ].copy()
    relocation_fixed["vs30_api_profile_mps"] = pd.to_numeric(
        relocation_fixed["vs30_api_profile_mps"], errors="coerce"
    )
    relocation_fixed["safe_vs30_proxy_mps"] = pd.to_numeric(
        relocation_fixed["safe_vs30_proxy_mps"], errors="coerce"
    )
    relocation_fixed = relocation_fixed.dropna(
        subset=["vs30_api_profile_mps", "safe_vs30_proxy_mps"]
    )
    if not relocation_fixed.empty:
        relocation_fixed["abs_delta"] = (
            relocation_fixed["vs30_api_profile_mps"]
            - relocation_fixed["safe_vs30_proxy_mps"]
        ).abs()
        relocation_fixed = relocation_fixed.sort_values(
            "abs_delta", ascending=False
        ).head(12)
        x = list(range(len(relocation_fixed)))
        fig, ax = plt.subplots(figsize=(8.4, 4.6))
        ax.bar(
            [i - 0.18 for i in x],
            relocation_fixed["vs30_api_profile_mps"],
            width=0.36,
            label="current station page",
            color="#b86f3f",
        )
        ax.bar(
            [i + 0.18 for i in x],
            relocation_fixed["safe_vs30_proxy_mps"],
            width=0.36,
            label="relocation-safe",
            color="#376e8a",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(relocation_fixed["station_code"], rotation=45, ha="right")
        ax.set_ylabel("VS30 proxy (m/s)")
        ax.set_title("Relocation Corrections to Candidate VS30 Proxies")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(figure_dir / "knet_relocation_vs30_correction.png", dpi=180)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    safe_knet, safe_profiles = assign_knet_safe_metadata(args.output_dir)
    write_figures(args.output_dir, safe_knet)
    write_report(args.output_dir, safe_knet, safe_profiles)


if __name__ == "__main__":
    main()
