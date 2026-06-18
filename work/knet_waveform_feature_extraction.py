#!/usr/bin/env python3
"""Stream K-NET accelerogram BSON and extract waveform-level features.

The accelerogram collection is large, so this script processes documents
sequentially and stores only compact per-component and per-record features.

Validation from the local data shows that, after demeaning, the maximum
absolute horizontal component matches header PGA_gal. Therefore the primary
time-domain features are computed from demeaned acceleration in gal-like units.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
import time
from pathlib import Path
from typing import Iterable

import bson
import numpy as np
import pandas as pd
from scipy.signal.windows import tukey


DEFAULT_KNET_DIR = Path("/Users/yojironoda/Downloads/s7rk7bj3zn-1/knet_1530")
DEFAULT_HEADER_INDEX = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_header_index.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)

FREQ_BANDS = [
    (0.2, 0.5),
    (0.5, 1.0),
    (1.0, 2.0),
    (2.0, 5.0),
    (5.0, 10.0),
    (10.0, 20.0),
]


def iter_bson_documents(path: Path) -> Iterable[dict]:
    with path.open("rb") as handle:
        while True:
            prefix = handle.read(4)
            if not prefix:
                break
            size = struct.unpack("<i", prefix)[0]
            payload = prefix + handle.read(size - 4)
            yield bson.loads(payload)


def safe_log10(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return float("nan")
    return float(np.log10(value))


def component_features(component_id: str, values: list[float], sampling_rate_hz: float) -> dict:
    record_id, component = component_id.rsplit(".", 1)
    arr = np.asarray(values, dtype=np.float64)
    npts = int(arr.size)
    if npts == 0:
        raise ValueError(f"Empty accelerogram for {component_id}")
    acc = arr - np.nanmean(arr)
    dt = 1.0 / float(sampling_rate_hz)
    duration_s = npts * dt
    abs_acc = np.abs(acc)

    # Velocity is a baseline-sensitive proxy. Demeaning makes it useful for
    # comparison, but publication-grade PGV should include robust correction.
    vel_proxy = np.cumsum(acc) * dt
    vel_proxy = vel_proxy - np.nanmean(vel_proxy)

    row = {
        "record_id": record_id,
        "component": component,
        "component_id": component_id,
        "npts": npts,
        "duration_s_waveform": duration_s,
        "sampling_rate_hz_waveform": sampling_rate_hz,
        "mean_before_demean": float(np.nanmean(arr)),
        "pga_component_gal": float(np.nanmax(abs_acc)),
        "rms_component_gal": float(np.sqrt(np.nanmean(acc**2))),
        "cav_component_proxy": float(np.nansum(abs_acc) * dt),
        "arias_component_proxy": float(np.nansum(acc**2) * dt),
        "pgv_component_proxy_cmps": float(np.nanmax(np.abs(vel_proxy))),
    }
    row["log10_pga_component_gal"] = safe_log10(row["pga_component_gal"])
    row["log10_rms_component_gal"] = safe_log10(row["rms_component_gal"])
    row["log10_cav_component_proxy"] = safe_log10(row["cav_component_proxy"])
    row["log10_arias_component_proxy"] = safe_log10(row["arias_component_proxy"])
    row["log10_pgv_component_proxy_cmps"] = safe_log10(row["pgv_component_proxy_cmps"])

    if npts >= 16:
        win = tukey(npts, alpha=0.08)
        spectrum = np.abs(np.fft.rfft(acc * win))
        freqs = np.fft.rfftfreq(npts, dt)
        for lo, hi in FREQ_BANDS:
            mask = (freqs >= lo) & (freqs < hi)
            amp = float(np.nanmedian(spectrum[mask])) if np.any(mask) else float("nan")
            key = f"spectral_amp_{lo:g}_{hi:g}_hz"
            row[key] = amp
            row[f"log10_{key}"] = safe_log10(amp)
    return row


def build_record_features(component_df: pd.DataFrame, header_df: pd.DataFrame) -> pd.DataFrame:
    pivot_cols = [
        "pga_component_gal",
        "rms_component_gal",
        "cav_component_proxy",
        "arias_component_proxy",
        "pgv_component_proxy_cmps",
    ]
    band_cols = [f"spectral_amp_{lo:g}_{hi:g}_hz" for lo, hi in FREQ_BANDS]
    pivot_cols += band_cols

    wide = component_df.pivot(index="record_id", columns="component", values=pivot_cols)
    wide.columns = [f"{metric}_{component}" for metric, component in wide.columns]
    wide = wide.reset_index()

    def cols_for(metric: str) -> tuple[str, str, str]:
        return f"{metric}_NS", f"{metric}_EW", f"{metric}_UD"

    ns, ew, ud = cols_for("pga_component_gal")
    wide["pga_horizontal_from_waveform_gal"] = wide[[ns, ew]].max(axis=1)
    wide["pga_3c_from_waveform_gal"] = wide[[ns, ew, ud]].max(axis=1)
    wide["pga_vertical_from_waveform_gal"] = wide[ud]

    for metric in [
        "rms_component_gal",
        "cav_component_proxy",
        "arias_component_proxy",
        "pgv_component_proxy_cmps",
    ]:
        ns, ew, ud = cols_for(metric)
        wide[f"{metric.replace('_component', '')}_horizontal"] = np.sqrt(
            (wide[ns] ** 2 + wide[ew] ** 2) / 2.0
        )
        wide[f"{metric.replace('_component', '')}_vertical"] = wide[ud]

    for lo, hi in FREQ_BANDS:
        metric = f"spectral_amp_{lo:g}_{hi:g}_hz"
        ns, ew, ud = cols_for(metric)
        wide[f"{metric}_horizontal"] = np.sqrt((wide[ns] ** 2 + wide[ew] ** 2) / 2.0)
        wide[f"{metric}_vertical"] = wide[ud]
        wide[f"log10_{metric}_horizontal"] = wide[f"{metric}_horizontal"].map(safe_log10)
        wide[f"log10_{metric}_vertical"] = wide[f"{metric}_vertical"].map(safe_log10)

    header_keep = header_df.copy()
    header_keep["_id"] = header_keep["_id"].astype(str)
    merged = header_keep.merge(wide, left_on="_id", right_on="record_id", how="left")
    merged["pga_header_waveform_ratio"] = (
        merged["pga_horizontal_from_waveform_gal"] / merged["PGA_gal"]
    )
    merged["log10_pga_horizontal_from_waveform_gal"] = merged[
        "pga_horizontal_from_waveform_gal"
    ].map(safe_log10)
    merged["log10_pga_3c_from_waveform_gal"] = merged["pga_3c_from_waveform_gal"].map(
        safe_log10
    )
    merged["log10_pga_vertical_from_waveform_gal"] = merged[
        "pga_vertical_from_waveform_gal"
    ].map(safe_log10)

    for col in [
        "rms_gal_horizontal",
        "cav_proxy_horizontal",
        "arias_proxy_horizontal",
        "pgv_proxy_cmps_horizontal",
        "rms_gal_vertical",
        "cav_proxy_vertical",
        "arias_proxy_vertical",
        "pgv_proxy_cmps_vertical",
    ]:
        if col in merged.columns:
            merged[f"log10_{col}"] = merged[col].map(safe_log10)
    return merged


def write_report(record_df: pd.DataFrame, component_df: pd.DataFrame, output_dir: Path) -> None:
    valid_ratio = record_df["pga_header_waveform_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
    lines = ["# K-NET Waveform Feature Extraction\n\n"]
    lines.append("## Inventory\n")
    lines.append(f"- Component feature rows: {len(component_df):,}\n")
    lines.append(f"- Record feature rows: {len(record_df):,}\n")
    lines.append(
        f"- Records with waveform features: {int(record_df['record_id'].notna().sum()):,} / {len(record_df):,}\n"
    )
    lines.append("\n## PGA consistency check\n")
    lines.append(
        "- Header PGA_gal should match the maximum absolute demeaned horizontal component.\n"
    )
    if len(valid_ratio):
        lines.append(f"- Median waveform/header PGA ratio: {valid_ratio.median():.8f}\n")
        lines.append(f"- 2.5%-97.5% range: {valid_ratio.quantile(0.025):.8f} to {valid_ratio.quantile(0.975):.8f}\n")
        lines.append(
            f"- Records with ratio within 1% of unity: {int(((valid_ratio > 0.99) & (valid_ratio < 1.01)).sum())} / {len(valid_ratio)}\n"
        )
    lines.append("\n## New waveform targets/features\n")
    lines.append("- Horizontal and vertical PGA from waveforms.\n")
    lines.append("- RMS acceleration, CAV proxy, Arias intensity proxy.\n")
    lines.append("- Baseline-sensitive PGV proxy from demeaned acceleration integration.\n")
    lines.append("- Fourier median amplitudes in 0.2-0.5, 0.5-1, 1-2, 2-5, 5-10, and 10-20 Hz bands.\n")
    lines.append("\n## Guardrails\n")
    lines.append("- PGV proxy is not final calibrated PGV; robust baseline correction is still required.\n")
    lines.append("- Spectral amplitudes are compact features, not response spectra.\n")
    lines.append("- Publication-grade spectral response should be added as a separate validated calculation.\n")
    (output_dir / "knet_waveform_feature_extraction.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--knet-dir", type=Path, default=DEFAULT_KNET_DIR)
    parser.add_argument("--header-index", type=Path, default=DEFAULT_HEADER_INDEX)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--max-components", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=5000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    header = pd.read_csv(args.header_index, dtype={"_id": str, "EventName": str})
    sr_by_record = header.set_index("_id")["SamplingRate_Hz"].to_dict()

    rows = []
    start = time.time()
    for idx, doc in enumerate(iter_bson_documents(args.knet_dir / "accelerogram.bson"), start=1):
        record_id, _ = doc["_id"].rsplit(".", 1)
        sampling_rate = float(sr_by_record.get(record_id, 100.0))
        rows.append(component_features(doc["_id"], doc["accelerogram"], sampling_rate))
        if args.progress_every and idx % args.progress_every == 0:
            elapsed = time.time() - start
            print(
                json.dumps(
                    {
                        "components": idx,
                        "elapsed_s": round(elapsed, 1),
                        "components_per_s": round(idx / max(elapsed, 1e-6), 1),
                    }
                ),
                flush=True,
            )
        if args.max_components and idx >= args.max_components:
            break

    component_df = pd.DataFrame(rows)
    component_df.to_csv(args.output_dir / "knet_component_waveform_features.csv", index=False)
    record_df = build_record_features(component_df, header)
    record_df.to_csv(args.output_dir / "knet_record_waveform_features.csv", index=False)
    write_report(record_df, component_df, args.output_dir)


if __name__ == "__main__":
    main()
