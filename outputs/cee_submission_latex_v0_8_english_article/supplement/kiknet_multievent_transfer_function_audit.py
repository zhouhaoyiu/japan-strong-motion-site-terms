#!/usr/bin/env python3
"""Multi-event KiK-net surface/downhole transfer-function audit.

The script accepts two local KiK-net layouts:

- converted miniSEED pairs under ``surface_data/waveforms`` and
  ``downhole_data/waveforms``;
- original NIED ASCII files with ``NS1/EW1/UD1`` and ``NS2/EW2/UD2`` channels.

For the original ASCII files the surface sensor is identified from the larger
station height in the paired headers, so the code does not assume that a fixed
suffix is always surface or downhole.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import read
from scipy.signal import detrend
from scipy.signal.windows import tukey


ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = ROOT / "outputs"
ARTICLE_DIR = OUTPUT_DIR / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT_DIR = ARTICLE_DIR / "supplement"
FIGURE_DIR = ARTICLE_DIR / "figures"

DEFAULT_MSEED_DIR = Path(
    os.environ.get("KIKNET_MSEED_DIR", str(ROOT / "work" / "external_data" / "kiknet_converted"))
)
DEFAULT_ASCII_DIR = Path(
    os.environ.get("KIKNET_ASCII_DIR", str(ROOT / "work" / "external_data" / "kiknet_ascii_event"))
)

OUT_ROWS = OUTPUT_DIR / "kiknet_multievent_transfer_functions.csv"
OUT_BINS = OUTPUT_DIR / "kiknet_multievent_transfer_function_bins.csv"
OUT_STATION = OUTPUT_DIR / "kiknet_multievent_transfer_function_station_summary.csv"
OUT_EVENT = OUTPUT_DIR / "kiknet_multievent_transfer_function_event_summary.csv"
OUT_AUDIT = OUTPUT_DIR / "kiknet_multievent_transfer_function_audit.md"
OUT_FIGURE = FIGURE_DIR / "figure16_kiknet_multievent_transfer_function.png"

COMPONENT_MAP = {"EW": "E", "NS": "N", "UD": "Z"}
FREQ_BANDS = [
    (0.2, 0.5),
    (0.5, 1.0),
    (1.0, 2.0),
    (2.0, 5.0),
    (5.0, 10.0),
    (10.0, 20.0),
]
FREQ_BIN_EDGES = np.geomspace(0.2, 20.0, 33)
ASCII_HEADER_KEYS = sorted(
    [
        "Origin Time",
        "Station Height(m)",
        "Sampling Freq(Hz)",
        "Duration Time(s)",
        "Scale Factor",
        "Max. Acc. (gal)",
        "Last Correction",
        "Station Code",
        "Station Lat.",
        "Station Long.",
        "Record Time",
        "Depth. (km)",
        "Lat.",
        "Long.",
        "Mag.",
        "Dir.",
        "Memo.",
    ],
    key=len,
    reverse=True,
)


@dataclass
class WaveformRecord:
    event_id: str
    station: str
    source_format: str
    surface: dict[str, np.ndarray]
    downhole: dict[str, np.ndarray]
    sampling_rate_hz: float
    surface_height_m: float | None = None
    downhole_height_m: float | None = None
    surface_path: str = ""
    downhole_path: str = ""


def safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0:
        return float("nan")
    return float(num / den)


def safe_log10(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return float("nan")
    return float(np.log10(value))


def preprocess(data: np.ndarray) -> np.ndarray:
    x = np.asarray(data, dtype=float)
    if len(x) < 2:
        return x
    x = detrend(x, type="linear")
    x = x - np.nanmean(x)
    return x


def amplitude_spectrum(data: np.ndarray, sampling_rate: float) -> tuple[np.ndarray, np.ndarray]:
    x = preprocess(data)
    if len(x) < 16:
        return np.array([]), np.array([])
    win = tukey(len(x), alpha=0.08)
    spec = np.abs(np.fft.rfft(x * win))
    freqs = np.fft.rfftfreq(len(x), d=1.0 / sampling_rate)
    return freqs, spec


def horizontal_series(comps: dict[str, np.ndarray]) -> np.ndarray | None:
    if "E" not in comps or "N" not in comps:
        return None
    n = min(len(comps["E"]), len(comps["N"]))
    return np.sqrt(preprocess(comps["E"][:n]) ** 2 + preprocess(comps["N"][:n]) ** 2)


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
    freqs: np.ndarray,
    surface_spec: np.ndarray,
    downhole_spec: np.ndarray,
    lo: float,
    hi: float,
) -> float:
    if len(freqs) == 0:
        return float("nan")
    mask = (freqs >= lo) & (freqs < hi)
    if not np.any(mask):
        return float("nan")
    eps = np.nanmedian(downhole_spec[mask]) * 1e-9 + 1e-12
    ratio = (surface_spec[mask] + eps) / (downhole_spec[mask] + eps)
    ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
    return float(np.nanmedian(ratio)) if len(ratio) else float("nan")


def binned_transfer_rows(
    event_id: str,
    station: str,
    source_format: str,
    freqs: np.ndarray,
    surface_spec: np.ndarray,
    downhole_spec: np.ndarray,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    if len(freqs) == 0:
        return rows
    global_slice = downhole_spec[(freqs >= 0.2) & (freqs <= 20.0)]
    eps_global = (float(np.nanmedian(global_slice)) if len(global_slice) else 0.0) * 1e-9 + 1e-12
    for lo, hi in zip(FREQ_BIN_EDGES[:-1], FREQ_BIN_EDGES[1:], strict=True):
        mask = (freqs >= lo) & (freqs < hi)
        if not np.any(mask):
            continue
        ratio = (surface_spec[mask] + eps_global) / (downhole_spec[mask] + eps_global)
        ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
        if len(ratio) == 0:
            continue
        rows.append(
            {
                "event_id": event_id,
                "station": station,
                "source_format": source_format,
                "freq_lo_hz": float(lo),
                "freq_hi_hz": float(hi),
                "freq_mid_hz": float(math.sqrt(lo * hi)),
                "transfer_ratio_median": float(np.nanmedian(ratio)),
                "log10_transfer_ratio_median": safe_log10(float(np.nanmedian(ratio))),
                "n_frequency_points": int(len(ratio)),
            }
        )
    return rows


def event_from_name(name: str) -> str:
    station_event = re.match(r"[A-Z]{3,4}H?\d{2}(\d{10,14})", name)
    if station_event:
        digits = station_event.group(1)
        return digits[:10]
    match = re.search(r"(\d{10,14})", name)
    return match.group(1) if match else "unknown"


def station_from_name(name: str) -> str:
    match = re.match(r"([A-Z]{3,4}H?\d{2})", name)
    if match:
        return match.group(1)
    return Path(name).stem[:6]


def stream_to_components(path: Path) -> tuple[dict[str, np.ndarray], float, int]:
    st = read(str(path))
    comps: dict[str, np.ndarray] = {}
    rates = set()
    npts = []
    for tr in st:
        tr = tr.copy()
        tr.detrend("linear")
        tr.detrend("demean")
        comp = tr.stats.channel[-1].upper()
        data = tr.data.astype(float)
        comps[comp] = data
        rates.add(float(tr.stats.sampling_rate))
        npts.append(int(tr.stats.npts))
    if len(rates) != 1:
        raise ValueError(f"Multiple sampling rates in {path}: {rates}")
    if len(set(npts)) != 1:
        min_len = min(npts)
        comps = {key: value[:min_len] for key, value in comps.items()}
    return comps, rates.pop(), min(len(value) for value in comps.values())


def mseed_records(kiknet_dir: Path) -> list[WaveformRecord]:
    surface_dir = kiknet_dir / "surface_data" / "waveforms"
    downhole_dir = kiknet_dir / "downhole_data" / "waveforms"
    if not surface_dir.exists() or not downhole_dir.exists():
        return []
    surface_files = {p.name: p for p in surface_dir.glob("*.mseed")}
    downhole_files = {p.name: p for p in downhole_dir.glob("*.mseed")}
    records: list[WaveformRecord] = []
    for name in sorted(set(surface_files) & set(downhole_files)):
        surface, sr_s, _ = stream_to_components(surface_files[name])
        downhole, sr_d, _ = stream_to_components(downhole_files[name])
        if abs(sr_s - sr_d) > 1e-6:
            raise ValueError(f"Sampling-rate mismatch in {name}: {sr_s} vs {sr_d}")
        records.append(
            WaveformRecord(
                event_id=event_from_name(name),
                station=station_from_name(name),
                source_format="converted_mseed",
                surface=surface,
                downhole=downhole,
                sampling_rate_hz=sr_s,
                surface_path=str(surface_files[name]),
                downhole_path=str(downhole_files[name]),
            )
        )
    return records


def parse_ascii_channel(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header: dict[str, str] = {}
    data_start = 0
    for i, line in enumerate(text):
        stripped = line.strip()
        if stripped.startswith("Memo."):
            data_start = i + 1
            break
        if not stripped:
            continue
        matched_key = None
        for key in ASCII_HEADER_KEYS:
            if stripped.startswith(key):
                matched_key = key
                header[key] = stripped[len(key) :].strip()
                break
        if matched_key is None:
            parts = re.split(r"\s{2,}", stripped, maxsplit=1)
            if len(parts) == 2:
                header[parts[0].strip()] = parts[1].strip()

    values: list[int] = []
    for line in text[data_start:]:
        values.extend(int(item) for item in re.findall(r"[-+]?\d+", line))
    data = np.asarray(values, dtype=float)

    scale_text = header.get("Scale Factor", "")
    scale_match = re.search(r"([0-9.]+)\(gal\)/([0-9.]+)", scale_text)
    if scale_match:
        data = data * (float(scale_match.group(1)) / float(scale_match.group(2)))

    freq_text = header.get("Sampling Freq(Hz)", "")
    freq_match = re.search(r"([0-9.]+)", freq_text)
    sampling_rate = float(freq_match.group(1)) if freq_match else float("nan")
    height = float(header.get("Station Height(m)", "nan"))
    suffix = path.suffix.upper().lstrip(".")
    comp = COMPONENT_MAP.get(suffix[:2])
    sensor = suffix[-1] if suffix else ""
    return {
        "path": path,
        "station": header.get("Station Code", station_from_name(path.name)),
        "event_id": event_from_name(path.name),
        "component": comp,
        "sensor": sensor,
        "sampling_rate_hz": sampling_rate,
        "height_m": height,
        "data": data,
    }


def ascii_records(ascii_dir: Path) -> list[WaveformRecord]:
    if not ascii_dir.exists():
        return []
    channel_files = [
        p
        for p in ascii_dir.iterdir()
        if p.is_file() and re.search(r"\.(NS|EW|UD)[12]$", p.name, flags=re.IGNORECASE)
    ]
    channels = [parse_ascii_channel(path) for path in channel_files]
    groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for channel in channels:
        groups.setdefault((str(channel["event_id"]), str(channel["station"])), []).append(channel)

    records: list[WaveformRecord] = []
    for (event_id, station), items in sorted(groups.items()):
        by_sensor: dict[str, dict[str, object]] = {}
        for item in items:
            sensor = str(item["sensor"])
            by_sensor.setdefault(
                sensor,
                {
                    "height_m": float(item["height_m"]),
                    "sampling_rate_hz": float(item["sampling_rate_hz"]),
                    "components": {},
                    "paths": [],
                },
            )
            comp = item["component"]
            if comp:
                by_sensor[sensor]["components"][comp] = item["data"]
            by_sensor[sensor]["paths"].append(str(item["path"]))
        complete = {
            sensor: values
            for sensor, values in by_sensor.items()
            if {"E", "N", "Z"}.issubset(set(values["components"]))
        }
        if len(complete) < 2:
            continue
        ordered = sorted(complete.items(), key=lambda pair: float(pair[1]["height_m"]))
        downhole_sensor, downhole_values = ordered[0]
        surface_sensor, surface_values = ordered[-1]
        sr_s = float(surface_values["sampling_rate_hz"])
        sr_d = float(downhole_values["sampling_rate_hz"])
        if abs(sr_s - sr_d) > 1e-6:
            raise ValueError(f"Sampling-rate mismatch for {event_id} {station}: {sr_s} vs {sr_d}")
        records.append(
            WaveformRecord(
                event_id=event_id,
                station=station,
                source_format="nied_ascii",
                surface={k: np.asarray(v, dtype=float) for k, v in surface_values["components"].items()},
                downhole={k: np.asarray(v, dtype=float) for k, v in downhole_values["components"].items()},
                sampling_rate_hz=sr_s,
                surface_height_m=float(surface_values["height_m"]),
                downhole_height_m=float(downhole_values["height_m"]),
                surface_path=";".join(surface_values["paths"]),
                downhole_path=";".join(downhole_values["paths"]),
            )
        )
    return records


def process_record(record: WaveformRecord) -> tuple[dict[str, object], list[dict[str, object]]]:
    s_h = horizontal_series(record.surface)
    d_h = horizontal_series(record.downhole)
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

    if "Z" in record.surface and "Z" in record.downhole:
        n = min(len(record.surface["Z"]), len(record.downhole["Z"]))
        z_surface = preprocess(record.surface["Z"][:n])
        z_downhole = preprocess(record.downhole["Z"][:n])
        vertical_pga_surface = float(np.nanmax(np.abs(z_surface)))
        vertical_pga_downhole = float(np.nanmax(np.abs(z_downhole)))
    else:
        vertical_pga_surface = vertical_pga_downhole = float("nan")

    freqs_s, spec_s = horizontal_spectrum(record.surface, record.sampling_rate_hz)
    freqs_d, spec_d = horizontal_spectrum(record.downhole, record.sampling_rate_hz)
    if len(freqs_s) != len(freqs_d):
        freqs = np.array([])
        spec_s = np.array([])
        spec_d = np.array([])
    else:
        freqs = freqs_s

    row: dict[str, object] = {
        "event_id": record.event_id,
        "station": record.station,
        "source_format": record.source_format,
        "sampling_rate_hz": record.sampling_rate_hz,
        "surface_npts": min(len(value) for value in record.surface.values()),
        "downhole_npts": min(len(value) for value in record.downhole.values()),
        "surface_height_m": record.surface_height_m,
        "downhole_height_m": record.downhole_height_m,
        "surface_horizontal_pga_proxy": horizontal_pga_surface,
        "downhole_horizontal_pga_proxy": horizontal_pga_downhole,
        "surface_horizontal_rms_proxy": horizontal_rms_surface,
        "downhole_horizontal_rms_proxy": horizontal_rms_downhole,
        "surface_vertical_pga_proxy": vertical_pga_surface,
        "downhole_vertical_pga_proxy": vertical_pga_downhole,
    }
    row["horizontal_pga_ratio"] = safe_ratio(horizontal_pga_surface, horizontal_pga_downhole)
    row["horizontal_rms_ratio"] = safe_ratio(horizontal_rms_surface, horizontal_rms_downhole)
    row["vertical_pga_ratio"] = safe_ratio(vertical_pga_surface, vertical_pga_downhole)
    row["log10_horizontal_pga_ratio"] = safe_log10(float(row["horizontal_pga_ratio"]))
    row["log10_horizontal_rms_ratio"] = safe_log10(float(row["horizontal_rms_ratio"]))
    row["log10_vertical_pga_ratio"] = safe_log10(float(row["vertical_pga_ratio"]))

    for lo, hi in FREQ_BANDS:
        ratio = median_band_ratio(freqs, spec_s, spec_d, lo, hi)
        key = f"spectral_ratio_{lo:g}_{hi:g}_hz"
        row[key] = ratio
        row[f"log10_{key}"] = safe_log10(ratio)

    bins = binned_transfer_rows(record.event_id, record.station, record.source_format, freqs, spec_s, spec_d)
    return row, bins


def summarize(rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    value_cols = [
        "horizontal_pga_ratio",
        "horizontal_rms_ratio",
        "vertical_pga_ratio",
        *[f"spectral_ratio_{lo:g}_{hi:g}_hz" for lo, hi in FREQ_BANDS],
    ]
    event_rows: list[dict[str, object]] = []
    for event_id, sub in rows.groupby("event_id", sort=True):
        item: dict[str, object] = {
            "event_id": event_id,
            "source_formats": ",".join(sorted(sub["source_format"].unique())),
            "n_station_events": int(len(sub)),
            "n_stations": int(sub["station"].nunique()),
        }
        for col in value_cols:
            vals = pd.to_numeric(sub[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            item[f"{col}_q25"] = float(vals.quantile(0.25)) if len(vals) else np.nan
            item[f"{col}_median"] = float(vals.quantile(0.50)) if len(vals) else np.nan
            item[f"{col}_q75"] = float(vals.quantile(0.75)) if len(vals) else np.nan
        event_rows.append(item)

    station_rows: list[dict[str, object]] = []
    for station, sub in rows.groupby("station", sort=True):
        item = {
            "station": station,
            "n_events": int(sub["event_id"].nunique()),
            "events": ",".join(sorted(sub["event_id"].astype(str).unique())),
        }
        for col in value_cols:
            vals = pd.to_numeric(sub[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            item[f"{col}_median"] = float(vals.median()) if len(vals) else np.nan
        station_rows.append(item)
    return pd.DataFrame(event_rows), pd.DataFrame(station_rows)


def save_figure(rows: pd.DataFrame, bins: pd.DataFrame, event_summary: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.9), constrained_layout=True)

    ax = axes[0, 0]
    labels = []
    data = []
    for event_id, sub in rows.groupby("event_id", sort=True):
        vals = pd.to_numeric(sub["log10_horizontal_pga_ratio"], errors="coerce").dropna()
        if len(vals):
            labels.append(event_id)
            data.append(vals.to_numpy())
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0.0, color="0.25", lw=0.9)
    ax.set_ylabel("log10 horizontal PGA ratio")
    ax.set_title("A Surface/downhole PGA ratio by event")
    ax.tick_params(axis="x", rotation=20)

    ax = axes[0, 1]
    band_cols = [f"log10_spectral_ratio_{lo:g}_{hi:g}_hz" for lo, hi in FREQ_BANDS]
    band_labels = [f"{lo:g}-{hi:g}" for lo, hi in FREQ_BANDS]
    data = [pd.to_numeric(rows[col], errors="coerce").dropna().to_numpy() for col in band_cols]
    ax.boxplot(data, tick_labels=band_labels, showfliers=False)
    ax.axhline(0.0, color="0.25", lw=0.9)
    ax.set_xlabel("Frequency band (Hz)")
    ax.set_ylabel("log10 transfer ratio")
    ax.set_title("B Frequency-band transfer ratios")
    ax.tick_params(axis="x", rotation=25)

    ax = axes[1, 0]
    colors = {"2512082315": "#4C78A8", "2601061018": "#C44E52"}
    for event_id, sub in bins.groupby("event_id", sort=True):
        curve = (
            sub.groupby("freq_mid_hz", as_index=False)["transfer_ratio_median"]
            .median()
            .sort_values("freq_mid_hz")
        )
        ax.plot(
            curve["freq_mid_hz"],
            curve["transfer_ratio_median"],
            lw=1.7,
            color=colors.get(str(event_id), None),
            label=f"{event_id} median",
        )
    ax.axhline(1.0, color="0.25", lw=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Surface/downhole transfer ratio")
    ax.set_title("C Median transfer functions")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    x = pd.to_numeric(rows["downhole_horizontal_pga_proxy"], errors="coerce")
    y = pd.to_numeric(rows["horizontal_pga_ratio"], errors="coerce")
    for event_id, sub in rows.assign(_x=x, _y=y).groupby("event_id", sort=True):
        mask = sub["_x"].gt(0) & sub["_y"].gt(0)
        ax.scatter(
            sub.loc[mask, "_x"],
            sub.loc[mask, "_y"],
            s=12,
            alpha=0.50,
            linewidths=0,
            color=colors.get(str(event_id), None),
            label=str(event_id),
        )
    ax.axhline(1.0, color="0.25", lw=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Downhole horizontal PGA proxy")
    ax.set_ylabel("Horizontal PGA ratio")
    ax.set_title("D Ratio versus input motion")
    ax.legend(frameon=False)

    fig.savefig(OUT_FIGURE, dpi=300)
    plt.close(fig)


def write_audit(rows: pd.DataFrame, bins: pd.DataFrame, event_summary: pd.DataFrame, station_summary: pd.DataFrame) -> None:
    finite = rows.replace([np.inf, -np.inf], np.nan)
    h = pd.to_numeric(finite["horizontal_pga_ratio"], errors="coerce").dropna()
    lines = [
        "# Multi-event KiK-net transfer-function audit",
        "",
        "## Data inventory",
        "",
        f"- Event-station pairs: {len(rows):,}.",
        f"- Events: {rows['event_id'].nunique():,}.",
        f"- Stations: {rows['station'].nunique():,}.",
        f"- Frequency-bin rows: {len(bins):,}.",
        "",
        "## Event coverage",
        "",
        "| event_id | source_formats | station-events | median horizontal PGA ratio | median 0.2-0.5 Hz ratio | median 1-2 Hz ratio |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for _, row in event_summary.iterrows():
        lines.append(
            f"| {row['event_id']} | {row['source_formats']} | {int(row['n_station_events'])} | "
            f"{row['horizontal_pga_ratio_median']:.3f} | "
            f"{row['spectral_ratio_0.2_0.5_hz_median']:.3f} | "
            f"{row['spectral_ratio_1_2_hz_median']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Combined transfer-function result",
            "",
            f"- Median horizontal PGA ratio: {h.median():.3f}.",
            f"- Interquartile horizontal PGA ratio: {h.quantile(0.25):.3f} to {h.quantile(0.75):.3f}.",
            f"- Event-station pairs with ratio > 1: {int((h > 1.0).sum())} / {len(h)}.",
            f"- Event-station pairs with ratio > 2: {int((h > 2.0).sum())} / {len(h)}.",
            "",
            "## Boundary",
            "",
            "This audit estimates observed surface/downhole transfer functions from two available KiK-net events. The 2026 ASCII event is converted with the NIED header scale factors. The 2025 converted miniSEED event is treated as a same-record amplitude-ratio data path, consistent with the existing local conversion. The result is a multi-event physical consistency check for station terms and frequency bands; it is not a nonlinear site-response model and does not replace a full KiK-net archive analysis.",
            "",
            "## Output tables",
            "",
            f"- `{OUT_ROWS.relative_to(ROOT)}`",
            f"- `{OUT_BINS.relative_to(ROOT)}`",
            f"- `{OUT_EVENT.relative_to(ROOT)}`",
            f"- `{OUT_STATION.relative_to(ROOT)}`",
        ]
    )
    OUT_AUDIT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_to_supplement(paths: list[Path]) -> None:
    SUPPLEMENT_DIR.mkdir(parents=True, exist_ok=True)
    for path in paths:
        shutil.copy2(path, SUPPLEMENT_DIR / path.name)


def main() -> None:
    global OUTPUT_DIR, OUT_ROWS, OUT_BINS, OUT_STATION, OUT_EVENT, OUT_AUDIT

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mseed-dir", type=Path, default=DEFAULT_MSEED_DIR)
    parser.add_argument("--ascii-dir", type=Path, default=DEFAULT_ASCII_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    OUTPUT_DIR = args.output_dir
    OUT_ROWS = OUTPUT_DIR / "kiknet_multievent_transfer_functions.csv"
    OUT_BINS = OUTPUT_DIR / "kiknet_multievent_transfer_function_bins.csv"
    OUT_STATION = OUTPUT_DIR / "kiknet_multievent_transfer_function_station_summary.csv"
    OUT_EVENT = OUTPUT_DIR / "kiknet_multievent_transfer_function_event_summary.csv"
    OUT_AUDIT = OUTPUT_DIR / "kiknet_multievent_transfer_function_audit.md"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    records = []
    records.extend(mseed_records(args.mseed_dir))
    records.extend(ascii_records(args.ascii_dir))
    if not records:
        raise RuntimeError("No KiK-net surface/downhole records found.")

    row_items: list[dict[str, object]] = []
    bin_items: list[dict[str, object]] = []
    for record in records:
        row, bins = process_record(record)
        row_items.append(row)
        bin_items.extend(bins)

    rows = pd.DataFrame(row_items).sort_values(["event_id", "station"]).reset_index(drop=True)
    bins = pd.DataFrame(bin_items).sort_values(["event_id", "station", "freq_mid_hz"]).reset_index(drop=True)
    event_summary, station_summary = summarize(rows)

    rows.to_csv(OUT_ROWS, index=False)
    bins.to_csv(OUT_BINS, index=False)
    event_summary.to_csv(OUT_EVENT, index=False)
    station_summary.to_csv(OUT_STATION, index=False)
    save_figure(rows, bins, event_summary)
    write_audit(rows, bins, event_summary, station_summary)
    copy_to_supplement([Path(__file__), OUT_ROWS, OUT_BINS, OUT_EVENT, OUT_STATION, OUT_AUDIT])

    print(
        "Multi-event KiK-net transfer-function audit complete: "
        f"{len(rows)} event-station pairs, {rows['event_id'].nunique()} events, "
        f"{rows['station'].nunique()} stations, {len(bins)} transfer-bin rows. "
        f"Median horizontal PGA ratio {rows['horizontal_pga_ratio'].median():.3f}."
    )


if __name__ == "__main__":
    main()
