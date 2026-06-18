#!/usr/bin/env python3
"""Summarize the official J-SHIS/NIED smrec_schema.tsv without materializing it.

The public flatfile subset stores a large strong-motion record table. This
script reads only selected columns in chunks and writes compact audit products
that can be cited in the project evidence chain.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"
DEFAULT_ZIP = PROJECT_ROOT / "work" / "external_data" / "jshis_gmf" / "flatfile_sub1-v2024.zip"

SMREC_USECOLS = [
    "smrec_id",
    "filebasename",
    "site_id",
    "eq_source_id",
    "length",
    "samplefreq",
    "maxaccrd050",
    "maxaccrd100",
    "maxvelrd050",
    "maxvelrd100",
    "maxaccv",
    "maxvelv",
    "sival",
    "sindo",
    "rsaccrd050d005t0030",
    "rsaccrd050d005t0100",
    "rsaccrd050d005t0300",
    "fault_dist",
    "lower_period",
    "upper_period",
    "multiple",
]

METRIC_COLUMNS = [
    "length",
    "samplefreq",
    "maxaccrd050",
    "maxaccrd100",
    "maxvelrd050",
    "maxvelrd100",
    "maxaccv",
    "maxvelv",
    "sival",
    "sindo",
    "rsaccrd050d005t0030",
    "rsaccrd050d005t0100",
    "rsaccrd050d005t0300",
    "fault_dist",
    "lower_period",
    "upper_period",
]

METRIC_LABELS = {
    "length": "Record length samples",
    "samplefreq": "Sample frequency Hz",
    "maxaccrd050": "Official maxacc rd050",
    "maxaccrd100": "Official maxacc rd100",
    "maxvelrd050": "Official maxvel rd050",
    "maxvelrd100": "Official maxvel rd100",
    "maxaccv": "Official vertical maxacc",
    "maxvelv": "Official vertical maxvel",
    "sival": "SI value",
    "sindo": "JMA intensity",
    "rsaccrd050d005t0030": "Official SA rd050 5pct T=0.3s",
    "rsaccrd050d005t0100": "Official SA rd050 5pct T=1.0s",
    "rsaccrd050d005t0300": "Official SA rd050 5pct T=3.0s",
    "fault_dist": "Fault distance km",
    "lower_period": "Lower valid period",
    "upper_period": "Upper valid period",
}

NETWORK_LABELS = {
    1: "K-NET",
    2: "KiK-net",
}


def read_small_schemas(zip_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with ZipFile(zip_path) as zf:
        with zf.open("site_schema.tsv") as handle:
            site = pd.read_csv(
                handle,
                sep="\t",
                usecols=["siteid2", "site_code", "site_name", "obs_network_id"],
            )
        with zf.open("source_schema.tsv") as handle:
            source = pd.read_csv(
                handle,
                sep="\t",
                usecols=[
                    "eq_source_id",
                    "jem_origin_time",
                    "jem_lat",
                    "jem_lon",
                    "jem_depth",
                    "mjma",
                    "mw",
                    "eq_event_name",
                ],
            )
    site["siteid2"] = pd.to_numeric(site["siteid2"], errors="coerce").astype("Int64")
    source["eq_source_id"] = pd.to_numeric(source["eq_source_id"], errors="coerce").astype("Int64")
    return site, source


def update_metric_summary(summary: dict[str, dict[str, float]], chunk: pd.DataFrame) -> None:
    for col in METRIC_COLUMNS:
        values = pd.to_numeric(chunk[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if values.empty:
            continue
        item = summary.setdefault(
            col,
            {"n": 0, "sum": 0.0, "min": np.inf, "max": -np.inf},
        )
        item["n"] += int(values.size)
        item["sum"] += float(values.sum())
        item["min"] = min(item["min"], float(values.min()))
        item["max"] = max(item["max"], float(values.max()))


def append_metric_sample(samples: list[pd.DataFrame], chunk: pd.DataFrame, rows_per_chunk: int) -> None:
    if len(chunk) <= rows_per_chunk:
        sampled = chunk[METRIC_COLUMNS].copy()
    else:
        sampled = chunk[METRIC_COLUMNS].sample(n=rows_per_chunk, random_state=20260605)
    samples.append(sampled)


def metric_summary_frame(summary: dict[str, dict[str, float]], sample_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in METRIC_COLUMNS:
        values = pd.to_numeric(sample_df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        item = summary.get(col, {"n": 0, "sum": np.nan, "min": np.nan, "max": np.nan})
        if item["n"]:
            mean = item["sum"] / item["n"]
        else:
            mean = np.nan
        quantiles = values.quantile([0.05, 0.5, 0.95]) if not values.empty else pd.Series(dtype=float)
        rows.append(
            {
                "metric": col,
                "label": METRIC_LABELS[col],
                "n_nonmissing": int(item["n"]),
                "min": item["min"],
                "mean": mean,
                "sample_q05": float(quantiles.get(0.05, np.nan)),
                "sample_median": float(quantiles.get(0.5, np.nan)),
                "sample_q95": float(quantiles.get(0.95, np.nan)),
                "max": item["max"],
                "quantile_note": "Quantiles estimated from deterministic per-chunk samples.",
            }
        )
    return pd.DataFrame(rows)


def counter_frame(counter: Counter, key_name: str, count_name: str) -> pd.DataFrame:
    return pd.DataFrame(
        [{key_name: key, count_name: value} for key, value in counter.items()],
        columns=[key_name, count_name],
    )


def summarize_smrec(zip_path: Path, output_dir: Path, chunksize: int, sample_rows_per_chunk: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    site, source = read_small_schemas(zip_path)
    site_lookup = site.drop_duplicates("siteid2").copy()
    source_lookup = source.drop_duplicates("eq_source_id").copy()

    metric_summary: dict[str, dict[str, float]] = {}
    metric_samples: list[pd.DataFrame] = []
    site_counts: Counter = Counter()
    source_counts: Counter = Counter()
    network_counts: Counter = Counter()
    samplefreq_counts: Counter = Counter()
    multiple_counts: Counter = Counter()
    total_rows = 0
    matched_site_rows = 0
    unmatched_site_ids: Counter = Counter()

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
                total_rows += len(chunk)
                chunk["site_id"] = pd.to_numeric(chunk["site_id"], errors="coerce").astype("Int64")
                chunk["siteid2"] = (chunk["site_id"] // 10).astype("Int64")
                chunk["eq_source_id"] = pd.to_numeric(chunk["eq_source_id"], errors="coerce").astype("Int64")

                joined = chunk.merge(site_lookup, on="siteid2", how="left")
                matched_site_rows += int(joined["site_code"].notna().sum())
                bad_ids = joined.loc[joined["site_code"].isna(), "site_id"].dropna().astype(int)
                unmatched_site_ids.update(bad_ids.tolist())

                site_counts.update(joined["siteid2"].dropna().astype(int).tolist())
                source_counts.update(joined["eq_source_id"].dropna().astype(int).tolist())
                network_counts.update(joined["obs_network_id"].dropna().astype(int).tolist())
                samplefreq_counts.update(pd.to_numeric(joined["samplefreq"], errors="coerce").dropna().tolist())
                multiple_counts.update(pd.to_numeric(joined["multiple"], errors="coerce").dropna().astype(int).tolist())
                update_metric_summary(metric_summary, joined)
                append_metric_sample(metric_samples, joined, sample_rows_per_chunk)

                print(
                    f"chunk={chunk_index} cumulative_rows={total_rows:,} "
                    f"matched_site_rows={matched_site_rows:,}",
                    flush=True,
                )

    sample_df = pd.concat(metric_samples, ignore_index=True) if metric_samples else pd.DataFrame(columns=METRIC_COLUMNS)
    metrics = metric_summary_frame(metric_summary, sample_df)
    metrics.to_csv(output_dir / "jshis_smrec_schema_sub1_metric_summary.csv", index=False)

    network_df = counter_frame(network_counts, "obs_network_id", "n_smrec_records")
    network_df["network_label"] = network_df["obs_network_id"].map(NETWORK_LABELS).fillna("unknown")
    network_df = network_df[["obs_network_id", "network_label", "n_smrec_records"]].sort_values("obs_network_id")
    network_df.to_csv(output_dir / "jshis_smrec_schema_sub1_network_counts.csv", index=False)

    samplefreq_df = counter_frame(samplefreq_counts, "samplefreq", "n_smrec_records").sort_values("samplefreq")
    samplefreq_df.to_csv(output_dir / "jshis_smrec_schema_sub1_samplefreq_counts.csv", index=False)

    multiple_df = counter_frame(multiple_counts, "multiple", "n_smrec_records").sort_values("multiple")
    multiple_df.to_csv(output_dir / "jshis_smrec_schema_sub1_multiple_counts.csv", index=False)

    top_sites = counter_frame(site_counts, "siteid2", "n_smrec_records").sort_values(
        "n_smrec_records", ascending=False
    )
    top_sites = top_sites.merge(site_lookup, on="siteid2", how="left")
    top_sites.head(100).to_csv(output_dir / "jshis_smrec_schema_sub1_top_sites.csv", index=False)

    top_events = counter_frame(source_counts, "eq_source_id", "n_smrec_records").sort_values(
        "n_smrec_records", ascending=False
    )
    top_events = top_events.merge(source_lookup, on="eq_source_id", how="left")
    top_events.head(100).to_csv(output_dir / "jshis_smrec_schema_sub1_top_events.csv", index=False)

    unmatched_df = counter_frame(unmatched_site_ids, "site_id", "n_smrec_records").sort_values(
        "n_smrec_records", ascending=False
    )
    unmatched_df.head(100).to_csv(output_dir / "jshis_smrec_schema_sub1_unmatched_site_ids.csv", index=False)

    write_report(
        output_dir=output_dir,
        zip_path=zip_path,
        total_rows=total_rows,
        matched_site_rows=matched_site_rows,
        site_counts=site_counts,
        source_counts=source_counts,
        network_df=network_df,
        samplefreq_df=samplefreq_df,
        multiple_df=multiple_df,
        metrics=metrics,
        top_sites=top_sites,
        top_events=top_events,
        unmatched_df=unmatched_df,
    )


def write_report(
    output_dir: Path,
    zip_path: Path,
    total_rows: int,
    matched_site_rows: int,
    site_counts: Counter,
    source_counts: Counter,
    network_df: pd.DataFrame,
    samplefreq_df: pd.DataFrame,
    multiple_df: pd.DataFrame,
    metrics: pd.DataFrame,
    top_sites: pd.DataFrame,
    top_events: pd.DataFrame,
    unmatched_df: pd.DataFrame,
) -> None:
    lines = ["# J-SHIS/NIED smrec_schema Subset Summary\n\n"]
    lines.append("## Archive\n")
    lines.append(f"- File: `{zip_path.name}`.\n")
    lines.append("- Member summarized: `smrec_schema.tsv`.\n")
    lines.append("- Parsing strategy: selected columns, chunked read, compact aggregate outputs only.\n\n")

    lines.append("## Record coverage\n")
    lines.append(f"- Strong-motion records: {total_rows:,}.\n")
    lines.append(f"- Unique `siteid2` values after `site_id // 10` mapping: {len(site_counts):,}.\n")
    lines.append(f"- Unique `eq_source_id` values: {len(source_counts):,}.\n")
    match_rate = matched_site_rows / total_rows if total_rows else np.nan
    lines.append(f"- Site-schema matched record rows: {matched_site_rows:,} ({match_rate:.3%}).\n")
    lines.append(f"- Unmatched `site_id` values observed: {unmatched_df['site_id'].nunique() if not unmatched_df.empty else 0:,}.\n\n")

    lines.append("## Network counts\n")
    for _, row in network_df.iterrows():
        lines.append(
            f"- obs_network_id={int(row['obs_network_id'])} ({row['network_label']}): "
            f"{int(row['n_smrec_records']):,} records.\n"
        )
    lines.append("\n## Sampling and flags\n")
    for _, row in samplefreq_df.iterrows():
        lines.append(f"- samplefreq={row['samplefreq']}: {int(row['n_smrec_records']):,} records.\n")
    for _, row in multiple_df.iterrows():
        lines.append(f"- multiple={int(row['multiple'])}: {int(row['n_smrec_records']):,} records.\n")
    lines.append("\n## Key metric ranges\n")
    key_metrics = [
        "maxaccrd050",
        "maxaccrd100",
        "fault_dist",
        "rsaccrd050d005t0030",
        "rsaccrd050d005t0100",
        "rsaccrd050d005t0300",
        "sindo",
    ]
    key = metrics[metrics["metric"].isin(key_metrics)].copy()
    for _, row in key.iterrows():
        lines.append(
            f"- {row['label']}: n={int(row['n_nonmissing']):,}, "
            f"min={row['min']:.4g}, mean={row['mean']:.4g}, "
            f"sample median={row['sample_median']:.4g}, max={row['max']:.4g}.\n"
        )
    lines.append("\n## Highest-record events\n")
    for _, row in top_events.head(10).iterrows():
        event_name = row["eq_event_name"] if pd.notna(row.get("eq_event_name")) else ""
        mjma = row["mjma"] if pd.notna(row.get("mjma")) else np.nan
        mw = row["mw"] if pd.notna(row.get("mw")) else np.nan
        lines.append(
            f"- eq_source_id={int(row['eq_source_id'])}: {int(row['n_smrec_records']):,} records, "
            f"JMA M={mjma}, Mw={mw}, time={row.get('jem_origin_time')}, name={event_name}.\n"
        )
    lines.append("\n## Highest-record sites\n")
    for _, row in top_sites.head(10).iterrows():
        lines.append(
            f"- {row.get('site_code')} (`siteid2`={int(row['siteid2'])}, "
            f"network={row.get('obs_network_id')}): {int(row['n_smrec_records']):,} records.\n"
        )
    lines.append("\n## Outputs\n")
    for name in [
        "jshis_smrec_schema_sub1_metric_summary.csv",
        "jshis_smrec_schema_sub1_network_counts.csv",
        "jshis_smrec_schema_sub1_samplefreq_counts.csv",
        "jshis_smrec_schema_sub1_multiple_counts.csv",
        "jshis_smrec_schema_sub1_top_events.csv",
        "jshis_smrec_schema_sub1_top_sites.csv",
        "jshis_smrec_schema_sub1_unmatched_site_ids.csv",
    ]:
        lines.append(f"- `{name}`\n")
    lines.append("\n## Guardrail\n")
    lines.append(
        "This is a compact audit of the public `flatfile_sub1-v2024.zip` archive, not a "
        "new model evaluation. The response-spectrum columns are reported by official "
        "J-SHIS field names; manuscript text should define the exact official measure "
        "before comparing them with locally computed K-NET features.\n"
    )

    (output_dir / "jshis_smrec_schema_sub1_summary.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunksize", type=int, default=200_000)
    parser.add_argument("--sample-rows-per-chunk", type=int, default=10_000)
    args = parser.parse_args()
    summarize_smrec(args.zip, args.output_dir, args.chunksize, args.sample_rows_per_chunk)


if __name__ == "__main__":
    main()
