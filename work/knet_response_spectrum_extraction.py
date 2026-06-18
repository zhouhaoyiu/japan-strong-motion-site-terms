#!/usr/bin/env python3
"""Extract 5%-damped response-spectrum features from K-NET accelerograms.

The script streams the K-NET accelerogram BSON file and computes pseudo
spectral acceleration (PSA) for each component at selected oscillator periods.
For a relative displacement oscillator

    u'' + 2*zeta*w*u' + w^2*u = -a_g

the pseudo spectral acceleration is w^2 * max(abs(u)).
"""

from __future__ import annotations

import argparse
import json
import struct
import time
from pathlib import Path
from typing import Iterable

import bson
import numpy as np
import pandas as pd
from scipy.signal import bilinear, lfilter


DEFAULT_KNET_DIR = Path("/Users/yojironoda/Downloads/s7rk7bj3zn-1/knet_1530")
DEFAULT_HEADER_INDEX = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs/knet_header_index.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/yojironoda/Documents/Codex/2026-06-05/eqtransformer/outputs"
)
DEFAULT_PERIODS = (0.3, 1.0, 3.0)


def iter_bson_documents(path: Path) -> Iterable[dict]:
    with path.open("rb") as handle:
        while True:
            prefix = handle.read(4)
            if not prefix:
                break
            size = struct.unpack("<i", prefix)[0]
            yield bson.loads(prefix + handle.read(size - 4))


def safe_log10(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return float("nan")
    return float(np.log10(value))


def response_spectrum_psa(
    acceleration: np.ndarray, sampling_rate_hz: float, periods: tuple[float, ...], damping: float
) -> dict[str, float]:
    """Compute pseudo spectral acceleration for one acceleration trace."""
    acc = np.asarray(acceleration, dtype=np.float64)
    acc = acc - np.nanmean(acc)
    fs = float(sampling_rate_hz)
    out: dict[str, float] = {}
    for period in periods:
        omega = 2.0 * np.pi / float(period)
        # Relative displacement transfer function from ground acceleration.
        b_cont = [-1.0]
        a_cont = [1.0, 2.0 * damping * omega, omega**2]
        b_dig, a_dig = bilinear(b_cont, a_cont, fs=fs)
        rel_disp = lfilter(b_dig, a_dig, acc)
        psa = float((omega**2) * np.nanmax(np.abs(rel_disp)))
        key = f"psa_{period:g}s_5pct"
        out[key] = psa
        out[f"log10_{key}"] = safe_log10(psa)
    return out


def component_response_features(
    component_id: str,
    values: list[float],
    sampling_rate_hz: float,
    periods: tuple[float, ...],
    damping: float,
) -> dict:
    record_id, component = component_id.rsplit(".", 1)
    row = {
        "record_id": record_id,
        "component": component,
        "component_id": component_id,
        "sampling_rate_hz": sampling_rate_hz,
        "npts": len(values),
    }
    row.update(response_spectrum_psa(np.asarray(values), sampling_rate_hz, periods, damping))
    return row


def build_record_response_features(component_df: pd.DataFrame, header_df: pd.DataFrame, periods: tuple[float, ...]) -> pd.DataFrame:
    metrics = [f"psa_{period:g}s_5pct" for period in periods]
    wide = component_df.pivot(index="record_id", columns="component", values=metrics)
    wide.columns = [f"{metric}_{component}" for metric, component in wide.columns]
    wide = wide.reset_index()

    for metric in metrics:
        ns = f"{metric}_NS"
        ew = f"{metric}_EW"
        ud = f"{metric}_UD"
        wide[f"{metric}_horizontal"] = wide[[ns, ew]].max(axis=1)
        wide[f"{metric}_3c"] = wide[[ns, ew, ud]].max(axis=1)
        wide[f"{metric}_vertical"] = wide[ud]
        wide[f"log10_{metric}_horizontal"] = wide[f"{metric}_horizontal"].map(safe_log10)
        wide[f"log10_{metric}_3c"] = wide[f"{metric}_3c"].map(safe_log10)
        wide[f"log10_{metric}_vertical"] = wide[f"{metric}_vertical"].map(safe_log10)

    header = header_df.copy()
    header["_id"] = header["_id"].astype(str)
    return header.merge(wide, left_on="_id", right_on="record_id", how="left")


def write_report(record_df: pd.DataFrame, component_df: pd.DataFrame, periods: tuple[float, ...], output_dir: Path) -> None:
    lines = ["# K-NET Response Spectrum Feature Extraction\n\n"]
    lines.append("## Inventory\n")
    lines.append(f"- Component response rows: {len(component_df):,}\n")
    lines.append(f"- Record response rows: {len(record_df):,}\n")
    lines.append(f"- Records with response features: {int(record_df['record_id'].notna().sum()):,} / {len(record_df):,}\n")
    lines.append("- Damping ratio: 5%\n")
    lines.append(f"- Periods: {', '.join(f'{p:g}s' for p in periods)}\n")
    lines.append("\n## Median horizontal PSA\n")
    for period in periods:
        col = f"psa_{period:g}s_5pct_horizontal"
        vals = record_df[col].replace([np.inf, -np.inf], np.nan).dropna()
        if len(vals):
            lines.append(f"- SA({period:g}s): median={vals.median():.4f}, 95th percentile={vals.quantile(0.95):.4f}\n")
    lines.append("\n## Guardrails\n")
    lines.append("- These are pseudo spectral accelerations from demeaned acceleration records.\n")
    lines.append("- Final publication should document baseline correction and compare against established strong-motion processing conventions.\n")
    (output_dir / "knet_response_spectrum_extraction.md").write_text(
        "".join(lines), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--knet-dir", type=Path, default=DEFAULT_KNET_DIR)
    parser.add_argument("--header-index", type=Path, default=DEFAULT_HEADER_INDEX)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--periods", type=float, nargs="+", default=list(DEFAULT_PERIODS))
    parser.add_argument("--damping", type=float, default=0.05)
    parser.add_argument("--max-components", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=5000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    periods = tuple(float(x) for x in args.periods)

    header = pd.read_csv(args.header_index, dtype={"_id": str, "EventName": str})
    sr_by_record = header.set_index("_id")["SamplingRate_Hz"].to_dict()

    rows = []
    start = time.time()
    for idx, doc in enumerate(iter_bson_documents(args.knet_dir / "accelerogram.bson"), start=1):
        record_id, _ = doc["_id"].rsplit(".", 1)
        sampling_rate = float(sr_by_record.get(record_id, 100.0))
        rows.append(
            component_response_features(
                doc["_id"], doc["accelerogram"], sampling_rate, periods, args.damping
            )
        )
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
    component_df.to_csv(args.output_dir / "knet_component_response_spectrum_features.csv", index=False)
    record_df = build_record_response_features(component_df, header, periods)
    record_df.to_csv(args.output_dir / "knet_record_response_spectrum_features.csv", index=False)
    write_report(record_df, component_df, periods, args.output_dir)


if __name__ == "__main__":
    main()
