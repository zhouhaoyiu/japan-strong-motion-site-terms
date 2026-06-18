#!/usr/bin/env python3
"""Integrate official J-SHIS/NIED site_schema.tsv with K-NET station terms."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats


DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)
DEFAULT_ZIP = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/work/external_data/jshis_gmf/flatfile_sub1-v2024.zip"
)
MATCH_THRESHOLD_KM = 0.2


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


def read_flatfile_schemas(zip_path: Path, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(zip_path) as zf:
        with zf.open("site_schema.tsv") as handle:
            site = pd.read_csv(handle, sep="\t")
        with zf.open("source_schema.tsv") as handle:
            source = pd.read_csv(handle, sep="\t")
    site.to_csv(output_dir / "jshis_site_schema_v2024_sub1.csv", index=False)
    source.to_csv(output_dir / "jshis_source_schema_v2024_sub1.csv", index=False)
    return site, source


def summarize_local_knet(output_dir: Path) -> pd.DataFrame:
    header = pd.read_csv(output_dir / "knet_header_index.csv")
    local = (
        header.groupby("StationCode")
        .agg(
            local_n_records=("StationCode", "size"),
            min_year=("Year", "min"),
            max_year=("Year", "max"),
            local_latitude=("StationLat", "median"),
            local_longitude=("StationLong", "median"),
            local_elevation_m=("StationHeight_m", "median"),
        )
        .reset_index()
        .rename(columns={"StationCode": "station_code"})
    )
    return local


def match_site_schema(site: pd.DataFrame, local: pd.DataFrame) -> pd.DataFrame:
    knet_site = site[site["obs_network_id"].eq(1)].copy()
    rows = []
    grouped = {code: rows.copy() for code, rows in knet_site.groupby("site_code")}
    for _, station in local.iterrows():
        candidates = grouped.get(station["station_code"])
        row = station.to_dict()
        row["jshis_match_found"] = False
        row["jshis_coordinate_delta_km"] = np.nan
        if candidates is None or candidates.empty:
            rows.append(row)
            continue
        candidate_rows = []
        for _, candidate in candidates.iterrows():
            distance = haversine_km(
                station["local_latitude"],
                station["local_longitude"],
                candidate["lat"],
                candidate["lon"],
            )
            candidate_rows.append((distance, candidate))
        distance, best = min(candidate_rows, key=lambda item: item[0])
        row["jshis_match_found"] = distance <= MATCH_THRESHOLD_KM
        row["jshis_coordinate_delta_km"] = distance
        for col in knet_site.columns:
            row[f"jshis_{col}"] = best[col]
        rows.append(row)
    matched = pd.DataFrame(rows)
    matched["jshis_match_status"] = np.select(
        [
            matched["jshis_match_found"],
            matched["jshis_site_code"].notna(),
        ],
        ["coordinate_matched", "station_code_only_coordinate_mismatch"],
        default="no_jshis_site_code",
    )
    return matched


def build_station_term_table(matched: pd.DataFrame, output_dir: Path, min_records: int) -> pd.DataFrame:
    adjusted = pd.read_csv(output_dir / "knet_adjusted_station_terms.csv").rename(
        columns={"n_records": "adjusted_term_n_records"}
    )
    matched_terms = matched.merge(
        adjusted,
        left_on="station_code",
        right_on="StationCode",
        how="inner",
    )
    matched_terms = matched_terms[
        matched_terms["adjusted_term_n_records"].ge(min_records)
        & matched_terms["jshis_match_found"]
    ].copy()
    for col in [
        "jshis_vs10",
        "jshis_vs20",
        "jshis_vs30",
        "jshis_avs30",
        "jshis_d1100",
        "jshis_d1400",
        "jshis_d1700",
        "jshis_d2100",
        "jshis_dbase",
        "regularized_term_log10_pga",
    ]:
        matched_terms[col] = pd.to_numeric(matched_terms[col], errors="coerce")
    for col in ["jshis_vs10", "jshis_vs20", "jshis_vs30", "jshis_avs30"]:
        matched_terms[f"log10_{col}"] = np.log10(matched_terms[col])
    matched_terms.to_csv(
        output_dir / f"jshis_knet_station_terms_site_schema_min{min_records}.csv",
        index=False,
    )
    return matched_terms


def official_multitarget_correlations(output_dir: Path, min_records: int) -> pd.DataFrame:
    site_terms = pd.read_csv(
        output_dir / f"jshis_knet_station_terms_site_schema_min{min_records}.csv"
    )
    multitarget = pd.read_csv(output_dir / "knet_multitarget_station_terms.csv")
    merged = multitarget.merge(
        site_terms[
            [
                "station_code",
                "jshis_vs10",
                "jshis_vs20",
                "jshis_vs30",
                "jshis_avs30",
                "jshis_d1100",
                "jshis_dbase",
                "jshis_coordinate_delta_km",
                "jshis_start_date",
                "jshis_end_date",
            ]
        ],
        left_on="StationCode",
        right_on="station_code",
        how="inner",
    )
    for col in ["jshis_vs10", "jshis_vs20", "jshis_vs30", "jshis_avs30", "jshis_d1100", "jshis_dbase"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")
        merged[f"log10_{col}"] = np.log10(merged[col])
    merged.to_csv(
        output_dir / f"jshis_knet_multitarget_station_terms_with_site_min{min_records}.csv",
        index=False,
    )

    features = [
        ("log10(VS10)", "log10_jshis_vs10"),
        ("log10(VS20)", "log10_jshis_vs20"),
        ("log10(VS30)", "log10_jshis_vs30"),
        ("log10(AVS30)", "log10_jshis_avs30"),
        ("log10(D1100)", "log10_jshis_d1100"),
        ("log10(Dbase)", "log10_jshis_dbase"),
    ]
    target_order = {
        label: idx
        for idx, label in enumerate(
            ["PGA H", "Arias H", "0.2-0.5 Hz", "0.5-1 Hz", "1-2 Hz", "2-5 Hz", "5-10 Hz", "10-20 Hz", "SA 0.3s", "SA 1.0s", "SA 3.0s"]
        )
    }
    rows = []
    for (target, label, group), sub in merged.groupby(["target", "target_label", "target_group"]):
        for feature_label, feature_col in features:
            data = sub[[feature_col, "regularized_station_term"]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(data) < 30:
                continue
            rho, pvalue = stats.spearmanr(data[feature_col], data["regularized_station_term"])
            rows.append(
                {
                    "target": target,
                    "target_label": label,
                    "target_group": group,
                    "target_order": target_order.get(label, 999),
                    "feature": feature_label,
                    "n": len(data),
                    "spearman_rho": float(rho),
                    "pvalue": float(pvalue),
                }
            )
    corr = pd.DataFrame(rows).sort_values(["target_order", "feature"])
    corr.to_csv(
        output_dir / f"jshis_knet_multitarget_station_site_correlations_min{min_records}.csv",
        index=False,
    )
    return corr


def write_figures(corr: pd.DataFrame, output_dir: Path, min_records: int) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    heat = corr.pivot(index="target_label", columns="feature", values="spearman_rho")
    label_order = ["PGA H", "Arias H", "0.2-0.5 Hz", "0.5-1 Hz", "1-2 Hz", "2-5 Hz", "5-10 Hz", "10-20 Hz", "SA 0.3s", "SA 1.0s", "SA 3.0s"]
    heat = heat.loc[[label for label in label_order if label in heat.index]]
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    im = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=25, ha="right")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            val = heat.iloc[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("Official J-SHIS Site Parameters vs K-NET Station Terms")
    fig.colorbar(im, ax=ax, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(
        figure_dir / f"jshis_knet_multitarget_station_site_correlation_heatmap_min{min_records}.png",
        dpi=180,
    )
    plt.close(fig)


def write_report(
    site: pd.DataFrame,
    source: pd.DataFrame,
    matched: pd.DataFrame,
    station_terms: pd.DataFrame,
    corr: pd.DataFrame,
    output_dir: Path,
    min_records: int,
) -> None:
    best = corr.reindex(corr["spearman_rho"].abs().sort_values(ascending=False).index).head(10)
    lines = ["# J-SHIS/NIED Flatfile Site Schema Integration\n\n"]
    lines.append("## Official files parsed\n")
    lines.append(f"- `site_schema.tsv`: {len(site):,} rows, {site['site_code'].nunique():,} unique site codes.\n")
    lines.append(f"- K-NET site rows in site_schema: {int(site['obs_network_id'].eq(1).sum()):,}.\n")
    lines.append(f"- `source_schema.tsv`: {len(source):,} rows.\n")
    lines.append("\n## K-NET coordinate-safe matching\n")
    lines.append(f"- Local K-NET stations checked: {len(matched):,}.\n")
    for status, count in matched["jshis_match_status"].value_counts().items():
        lines.append(f"- {status}: {int(count):,}.\n")
    lines.append(
        f"- Station terms with n_records >= {min_records} and official coordinate-matched site rows: {len(station_terms):,}.\n"
    )
    lines.append("\n## Strongest official site-parameter associations\n")
    for _, row in best.iterrows():
        lines.append(
            f"- {row['target_label']} vs {row['feature']}: rho={row['spearman_rho']:.3f}, n={int(row['n'])}.\n"
        )
    lines.append("\n## Interpretation\n")
    lines.append(
        "The official J-SHIS site_schema confirms the earlier Kyoshin-page result: "
        "site-parameter associations are weak for scalar PGA and stronger for "
        "frequency-dependent Fourier/SA station terms. AVS30 and basin-depth fields "
        "are now available for manuscript-grade interpretation, while many K-NET "
        "rows still have missing borehole-derived VS30.\n\n"
    )
    lines.append("## Outputs\n")
    lines.append("- `jshis_site_schema_v2024_sub1.csv`\n")
    lines.append("- `jshis_source_schema_v2024_sub1.csv`\n")
    lines.append(f"- `jshis_knet_site_schema_match.csv`\n")
    lines.append(f"- `jshis_knet_station_terms_site_schema_min{min_records}.csv`\n")
    lines.append(f"- `jshis_knet_multitarget_station_site_correlations_min{min_records}.csv`\n")
    lines.append(
        f"- `figures/jshis_knet_multitarget_station_site_correlation_heatmap_min{min_records}.png`\n\n"
    )
    lines.append("## Guardrail\n")
    lines.append(
        "The downloaded file is `flatfile_sub1-v2024.zip`, a public M>=5 / fault-distance-limited subset. "
        "Its site_schema appears to include the full station metadata, but final claims should cite the "
        "official J-SHIS/NIED flatfile version and state exactly which subset archive was used.\n"
    )
    (output_dir / f"jshis_site_schema_integration_min{min_records}.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-records", type=int, default=10)
    args = parser.parse_args()

    site, source = read_flatfile_schemas(args.zip, args.output_dir)
    local = summarize_local_knet(args.output_dir)
    matched = match_site_schema(site, local)
    matched.to_csv(args.output_dir / "jshis_knet_site_schema_match.csv", index=False)
    station_terms = build_station_term_table(matched, args.output_dir, args.min_records)
    corr = official_multitarget_correlations(args.output_dir, args.min_records)
    write_figures(corr, args.output_dir, args.min_records)
    write_report(site, source, matched, station_terms, corr, args.output_dir, args.min_records)


if __name__ == "__main__":
    main()
