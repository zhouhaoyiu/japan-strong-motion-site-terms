#!/usr/bin/env python3
"""KiK-net surface/downhole amplification proxy audit.

The local converted KiK-net data contain one event with paired surface and
downhole miniSEED records. This script estimates station-level surface/downhole
amplification proxies from waveform amplitudes and spectra.

Important scope:
- Waveforms appear to be in count-like units with offsets, so traces are
  demeaned/detrended before ratios are computed.
- Ratios are meaningful as same-record surface/downhole proxies, but not as
  calibrated physical accelerations unless the conversion metadata are added.
- A single event cannot establish nonlinear site response. It can identify
  candidate stations and frequency bands for full KiK-net follow-up.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import read
from scipy.signal.windows import tukey


DEFAULT_KIKNET_DIR = Path("/Users/yojironoda/Desktop/iemeq/kiknet_converted")
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


def station_from_filename(path: Path) -> str:
    return path.name[:6]


def common_mseed_files(kiknet_dir: Path) -> list[tuple[str, Path, Path]]:
    surface = kiknet_dir / "surface_data" / "waveforms"
    downhole = kiknet_dir / "downhole_data" / "waveforms"
    surface_files = {p.name: p for p in surface.glob("*.mseed")}
    downhole_files = {p.name: p for p in downhole.glob("*.mseed")}
    names = sorted(set(surface_files) & set(downhole_files))
    return [(station_from_filename(Path(name)), surface_files[name], downhole_files[name]) for name in names]


def stream_to_components(path: Path) -> tuple[dict[str, np.ndarray], float, int]:
    st = read(str(path))
    comps: dict[str, np.ndarray] = {}
    sampling_rates = set()
    npts = []
    for tr in st:
        tr = tr.copy()
        tr.detrend("linear")
        tr.detrend("demean")
        comp = tr.stats.channel[-1].upper()
        data = tr.data.astype(float)
        comps[comp] = data
        sampling_rates.add(float(tr.stats.sampling_rate))
        npts.append(int(tr.stats.npts))
    if len(sampling_rates) != 1:
        raise ValueError(f"Multiple sampling rates in {path}: {sampling_rates}")
    if len(set(npts)) != 1:
        min_len = min(npts)
        comps = {k: v[:min_len] for k, v in comps.items()}
        npts = [min_len]
    return comps, sampling_rates.pop(), npts[0]


def horizontal_series(comps: dict[str, np.ndarray]) -> np.ndarray | None:
    if "E" not in comps or "N" not in comps:
        return None
    n = min(len(comps["E"]), len(comps["N"]))
    return np.sqrt(comps["E"][:n] ** 2 + comps["N"][:n] ** 2)


def amplitude_spectrum(data: np.ndarray, sampling_rate: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(data, dtype=float)
    if len(x) < 16:
        return np.array([]), np.array([])
    x = x - np.nanmean(x)
    win = tukey(len(x), alpha=0.08)
    spec = np.abs(np.fft.rfft(x * win))
    freqs = np.fft.rfftfreq(len(x), d=1.0 / sampling_rate)
    return freqs, spec


def horizontal_spectrum(comps: dict[str, np.ndarray], sampling_rate: float) -> tuple[np.ndarray, np.ndarray]:
    if "E" not in comps or "N" not in comps:
        return np.array([]), np.array([])
    n = min(len(comps["E"]), len(comps["N"]))
    f_e, s_e = amplitude_spectrum(comps["E"][:n], sampling_rate)
    f_n, s_n = amplitude_spectrum(comps["N"][:n], sampling_rate)
    if len(f_e) == 0 or len(f_n) == 0 or len(f_e) != len(f_n):
        return np.array([]), np.array([])
    return f_e, np.sqrt((s_e**2 + s_n**2) / 2.0)


def median_band_ratio(
    freqs: np.ndarray, surface_spec: np.ndarray, downhole_spec: np.ndarray, lo: float, hi: float
) -> float:
    if len(freqs) == 0:
        return float("nan")
    mask = (freqs >= lo) & (freqs < hi)
    if not np.any(mask):
        return float("nan")
    eps = np.nanmedian(downhole_spec[mask]) * 1e-9 + 1e-12
    ratio = (surface_spec[mask] + eps) / (downhole_spec[mask] + eps)
    ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
    if len(ratio) == 0:
        return float("nan")
    return float(np.nanmedian(ratio))


def process_pair(station: str, surface_path: Path, downhole_path: Path) -> dict:
    s_comps, s_sr, s_npts = stream_to_components(surface_path)
    d_comps, d_sr, d_npts = stream_to_components(downhole_path)
    if abs(s_sr - d_sr) > 1e-6:
        raise ValueError(f"Sampling-rate mismatch for {station}: {s_sr} vs {d_sr}")

    s_h = horizontal_series(s_comps)
    d_h = horizontal_series(d_comps)
    if s_h is None or d_h is None:
        horizontal_pga_surface = horizontal_pga_downhole = float("nan")
        horizontal_rms_surface = horizontal_rms_downhole = float("nan")
    else:
        n = min(len(s_h), len(d_h))
        s_h = s_h[:n]
        d_h = d_h[:n]
        horizontal_pga_surface = float(np.nanmax(np.abs(s_h)))
        horizontal_pga_downhole = float(np.nanmax(np.abs(d_h)))
        horizontal_rms_surface = float(np.sqrt(np.nanmean(s_h**2)))
        horizontal_rms_downhole = float(np.sqrt(np.nanmean(d_h**2)))

    if "Z" in s_comps and "Z" in d_comps:
        n = min(len(s_comps["Z"]), len(d_comps["Z"]))
        z_pga_surface = float(np.nanmax(np.abs(s_comps["Z"][:n])))
        z_pga_downhole = float(np.nanmax(np.abs(d_comps["Z"][:n])))
    else:
        z_pga_surface = z_pga_downhole = float("nan")

    freqs_s, spec_s = horizontal_spectrum(s_comps, s_sr)
    freqs_d, spec_d = horizontal_spectrum(d_comps, d_sr)
    if len(freqs_s) != len(freqs_d):
        freqs = np.array([])
        spec_s = np.array([])
        spec_d = np.array([])
    else:
        freqs = freqs_s

    row = {
        "station": station,
        "file": surface_path.name,
        "sampling_rate_hz": s_sr,
        "surface_npts": s_npts,
        "downhole_npts": d_npts,
        "surface_horizontal_pga_proxy": horizontal_pga_surface,
        "downhole_horizontal_pga_proxy": horizontal_pga_downhole,
        "surface_horizontal_rms_proxy": horizontal_rms_surface,
        "downhole_horizontal_rms_proxy": horizontal_rms_downhole,
        "surface_vertical_pga_proxy": z_pga_surface,
        "downhole_vertical_pga_proxy": z_pga_downhole,
    }
    row["horizontal_pga_ratio"] = safe_ratio(
        row["surface_horizontal_pga_proxy"], row["downhole_horizontal_pga_proxy"]
    )
    row["horizontal_rms_ratio"] = safe_ratio(
        row["surface_horizontal_rms_proxy"], row["downhole_horizontal_rms_proxy"]
    )
    row["vertical_pga_ratio"] = safe_ratio(
        row["surface_vertical_pga_proxy"], row["downhole_vertical_pga_proxy"]
    )
    row["log10_horizontal_pga_ratio"] = safe_log10(row["horizontal_pga_ratio"])
    row["log10_horizontal_rms_ratio"] = safe_log10(row["horizontal_rms_ratio"])
    row["log10_vertical_pga_ratio"] = safe_log10(row["vertical_pga_ratio"])

    for lo, hi in FREQ_BANDS:
        ratio = median_band_ratio(freqs, spec_s, spec_d, lo, hi)
        key = f"spectral_ratio_{lo:g}_{hi:g}_hz"
        row[key] = ratio
        row[f"log10_{key}"] = safe_log10(ratio)

    return row


def safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0:
        return float("nan")
    return float(num / den)


def safe_log10(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return float("nan")
    return float(np.log10(value))


def save_figures(df: pd.DataFrame, output_dir: Path) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    valid = df.replace([np.inf, -np.inf], np.nan)

    plt.figure(figsize=(7.2, 4.4))
    x = valid["log10_horizontal_pga_ratio"].dropna()
    plt.hist(x, bins=35, color="#4c78a8")
    plt.axvline(0, color="k", linewidth=1)
    plt.xlabel("log10(surface/downhole horizontal PGA proxy)")
    plt.ylabel("Station count")
    plt.title("KiK-net surface/downhole horizontal PGA ratio")
    plt.tight_layout()
    plt.savefig(figures / "kiknet_horizontal_pga_ratio_histogram.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.8))
    x = valid["downhole_horizontal_pga_proxy"]
    y = valid["horizontal_pga_ratio"]
    mask = x.notna() & y.notna() & (x > 0) & (y > 0)
    plt.scatter(x[mask], y[mask], s=24, alpha=0.72, color="#59a14f", edgecolors="k", linewidths=0.2)
    plt.xscale("log")
    plt.yscale("log")
    plt.axhline(1, color="k", linewidth=1)
    plt.xlabel("Downhole horizontal PGA proxy (detrended counts)")
    plt.ylabel("Surface/downhole horizontal PGA ratio")
    plt.title("Amplification proxy versus downhole motion")
    plt.tight_layout()
    plt.savefig(figures / "kiknet_ratio_vs_downhole_motion.png", dpi=220)
    plt.close()

    band_cols = [f"log10_spectral_ratio_{lo:g}_{hi:g}_hz" for lo, hi in FREQ_BANDS]
    labels = [f"{lo:g}-{hi:g}" for lo, hi in FREQ_BANDS]
    data = [valid[col].dropna().to_numpy() for col in band_cols]
    plt.figure(figsize=(8.2, 4.8))
    plt.boxplot(data, tick_labels=labels, showfliers=False)
    plt.axhline(0, color="k", linewidth=1)
    plt.xlabel("Frequency band (Hz)")
    plt.ylabel("log10(surface/downhole spectral ratio)")
    plt.title("Frequency-dependent KiK-net amplification proxy")
    plt.tight_layout()
    plt.savefig(figures / "kiknet_spectral_ratio_bands.png", dpi=220)
    plt.close()

    top = valid.sort_values("log10_horizontal_pga_ratio", ascending=False).head(20)
    plt.figure(figsize=(8.5, 5.2))
    plt.bar(top["station"], top["log10_horizontal_pga_ratio"], color="#f28e2b")
    plt.xticks(rotation=70, ha="right")
    plt.ylabel("log10(horizontal PGA ratio)")
    plt.title("Top KiK-net surface/downhole PGA ratios")
    plt.tight_layout()
    plt.savefig(figures / "kiknet_top_pga_ratio_stations.png", dpi=220)
    plt.close()


def write_report(df: pd.DataFrame, output_dir: Path) -> None:
    valid = df.replace([np.inf, -np.inf], np.nan)
    ratio = valid["horizontal_pga_ratio"].dropna()
    log_ratio = valid["log10_horizontal_pga_ratio"].dropna()
    rms_ratio = valid["horizontal_rms_ratio"].dropna()
    lines = ["# KiK-net Surface/Downhole Amplification Proxy Audit\n\n"]
    lines.append("## Data inventory\n")
    lines.append(f"- Paired surface/downhole miniSEED records: {len(df):,}\n")
    lines.append(f"- Sampling rates: {sorted(valid['sampling_rate_hz'].dropna().unique().tolist())}\n")
    lines.append(
        f"- Median record length: {int(valid['surface_npts'].median()) if len(valid) else 0:,} samples\n"
    )
    lines.append("\n## Horizontal amplification proxy\n")
    if len(ratio):
        lines.append(f"- Median surface/downhole horizontal PGA ratio: {ratio.median():.3f}\n")
        lines.append(f"- Interquartile range: {ratio.quantile(0.25):.3f} to {ratio.quantile(0.75):.3f}\n")
        lines.append(f"- Stations with ratio > 1: {int((ratio > 1).sum())} / {len(ratio)}\n")
        lines.append(f"- Stations with ratio > 2: {int((ratio > 2).sum())} / {len(ratio)}\n")
        lines.append(f"- Median horizontal RMS ratio: {rms_ratio.median():.3f}\n")
        if len(log_ratio) > 2:
            corr = valid[["log10_horizontal_pga_ratio", "downhole_horizontal_pga_proxy"]].dropna()
            corr = corr[corr["downhole_horizontal_pga_proxy"] > 0].copy()
            corr["log10_downhole_horizontal_pga_proxy"] = np.log10(
                corr["downhole_horizontal_pga_proxy"]
            )
            rho = corr["log10_horizontal_pga_ratio"].corr(
                corr["log10_downhole_horizontal_pga_proxy"], method="spearman"
            )
            lines.append(
                f"- Spearman correlation with downhole-motion level: {rho:.3f} "
                "(single-event diagnostic only)\n"
            )

    lines.append("\n## Frequency-dependent ratios\n")
    for lo, hi in FREQ_BANDS:
        col = f"spectral_ratio_{lo:g}_{hi:g}_hz"
        vals = valid[col].dropna()
        if len(vals):
            lines.append(f"- {lo:g}-{hi:g} Hz median spectral ratio: {vals.median():.3f}\n")

    lines.append("\n## Highest horizontal PGA ratios\n")
    cols = [
        "station",
        "horizontal_pga_ratio",
        "horizontal_rms_ratio",
        "spectral_ratio_1_2_hz",
        "spectral_ratio_2_5_hz",
        "spectral_ratio_5_10_hz",
    ]
    for _, row in valid.sort_values("horizontal_pga_ratio", ascending=False).head(12).iterrows():
        parts = [f"- {row['station']}: PGA ratio={row['horizontal_pga_ratio']:.2f}"]
        for col in cols[2:]:
            if col in row and np.isfinite(row[col]):
                parts.append(f"{col.replace('spectral_ratio_', 'SR ').replace('_hz', 'Hz')}={row[col]:.2f}")
        lines.append(", ".join(parts) + "\n")

    lines.append("\n## Interpretation guardrails\n")
    lines.append(
        "- These are surface/downhole amplification proxies from one event and count-like waveforms after detrending.\n"
    )
    lines.append(
        "- They support candidate-site discovery and frequency-band targeting, not a final nonlinear site-response claim.\n"
    )
    lines.append(
        "- The next publication-grade step is to merge full KiK-net multi-event data with station metadata such as borehole depth and Vs profiles.\n"
    )
    (output_dir / "kiknet_surface_downhole_audit.md").write_text("".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kiknet-dir", type=Path, default=DEFAULT_KIKNET_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for station, surface_path, downhole_path in common_mseed_files(args.kiknet_dir):
        rows.append(process_pair(station, surface_path, downhole_path))
    df = pd.DataFrame(rows).sort_values("station")
    df.to_csv(args.output_dir / "kiknet_surface_downhole_ratios.csv", index=False)
    save_figures(df, args.output_dir)
    write_report(df, args.output_dir)


if __name__ == "__main__":
    main()
