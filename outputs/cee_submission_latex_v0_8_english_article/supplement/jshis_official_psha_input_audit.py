#!/usr/bin/env python3
"""Audit public J-SHIS PSHA parameter inputs for reproducibility boundaries.

The audit records what is publicly available in the official J-SHIS PSHM
parameter packages and what remains outside a direct official full-PSHA rerun.
It does not redistribute raw J-SHIS zip files and does not claim to reproduce
the official national hazard maps.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import shutil
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work" / "jshis_full_psha_audit"
DOWNLOAD_DIR = WORK / "downloaded"
EXTRACT_DIR = WORK / "extracted"
OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = OUTPUT_DIR / "figures"

VERSION = "Y2024"
CSV_ZIP = f"P-{VERSION}-PRM.zip"
SHAPE_ZIP = f"P-{VERSION}-PRM-SHAPE.zip"
CSV_URL = f"https://www.j-shis.bosai.go.jp/map/JSHIS2/data/P/{VERSION}/{CSV_ZIP}"
SHAPE_URL = f"https://www.j-shis.bosai.go.jp/map/JSHIS2/data/P/{VERSION}/PRM/{SHAPE_ZIP}"
DOWNLOAD_PAGE_URL = "https://www.j-shis.bosai.go.jp/map/JSHIS2/download.html?lang=en"
API_LIST_URL = "https://www.j-shis.bosai.go.jp/en/api-list"
DATA_RULE_URL = "https://www.j-shis.bosai.go.jp/map/JSHIS2/data/DOC/DataFileRule/A-RULES_en.pdf"


@dataclass(frozen=True)
class Category:
    key: str
    label: str
    role: str
    fields: str
    status: str
    paper_use: str


CATEGORIES = {
    "activity_characteristic": Category(
        "activity_characteristic",
        "Seismic activity parameters for specified faults",
        "Recurrence process, mean recurrence interval, latest-event time, variance, and 30/50-year occurrence probabilities for specified active-fault and subduction sources.",
        "CODE, PROC, AVRACT, NEWACT, ALPHA, P_T30, P_T50, NAME",
        "Public input. Needs conversion from J-SHIS CSV conventions to an executable PSHA source model.",
        "Supports a documented official-source input-chain audit.",
    ),
    "activity_background": Category(
        "activity_background",
        "Occurrence-frequency grids for sources without specified faults",
        "Mesh-based occurrence frequencies, b-values, minimum magnitudes, zone identifiers, depths, strikes, and dips for background or area sources.",
        "MNO, JLG, JLA, WLG, WLA, FRQ, BVL, MMN, ANO, DEP, STR, DIP",
        "Public input. Needs gridded-source translation and validation against J-SHIS calculation settings.",
        "Documents the public background-source input layer.",
    ),
    "source_geometry": Category(
        "source_geometry",
        "Fault/source geometry",
        "Rectangular, non-rectangular, and discretized source geometries for specified and area sources.",
        "Source code, magnitude or frequency block, point/fault-plane coordinates, depth, length, width, strike, dip",
        "Public input in CSV and shapefile forms. Geometry still needs source-model conversion for a full rerun.",
        "Shows that official source geometry is downloadable and versioned.",
    ),
    "zoning_area": Category(
        "zoning_area",
        "Zoning-area polygons",
        "Area boundaries used by source-frequency models.",
        "JLON, JLAT, WLON, WLAT",
        "Public input. Needed to reconstruct area-source occurrence models.",
        "Documents spatial zoning used by background-source files.",
    ),
    "inter_intra_ratio": Category(
        "inter_intra_ratio",
        "Interplate/intraplate occurrence ratio",
        "Ratios used to split plate-interface and intraslab occurrence components.",
        "EQCODE, ANO, INTERR, INTRAR",
        "Public input. Must be mapped into the executable source logic.",
        "Documents plate-source component ratios.",
    ),
    "attenuation_formula": Category(
        "attenuation_formula",
        "Attenuation-relation parameter codes",
        "J-SHIS code mapping for earthquake type, source-position type, magnitude-conversion type, and crustal-type flag.",
        "EQCODE, EQTYPE, SPTYPE, MTTYPE, CRTYPE",
        "Public input. The numerical GMPE implementation and calculation switches must be verified before claiming a full official rerun.",
        "Defines the remaining bridge between official source inputs and ground-motion calculation.",
    ),
    "plate_shape": Category(
        "plate_shape",
        "Pacific/Philippine Sea plate depth surfaces",
        "Gridded plate-surface depth points for plate-interface and intraslab source treatment.",
        "MNO, JLG, JLA, WLG, WLA, DEP",
        "Public input. Needs validated interpolation and source-placement logic.",
        "Documents public plate-geometry inputs.",
    ),
    "eqthr": Category(
        "eqthr",
        "EQTHR activity and geometry",
        "Parameters and geometry for earthquakes whose traces are hard to recognize from surface evidence.",
        "Fault code, recurrence interval, magnitude bounds, b-value, and fault geometry",
        "Public input. Needs conversion and validation.",
        "Documents the active-fault supplement layer.",
    ),
    "other": Category(
        "other",
        "Other J-SHIS PSHM parameter files",
        "Files outside the main categories detected by the audit script.",
        "Varies by file",
        "Review before using in a calculation.",
        "Inventory only.",
    ),
}


def category_for(filename: str) -> str:
    stem = filename.replace(f"P-{VERSION}-PRM-", "").replace(".csv", "")
    if stem.startswith("ACT_AVR_") or stem.startswith("ACT_MAX_"):
        return "activity_characteristic"
    if stem.endswith("EQTHR"):
        return "eqthr"
    if stem.startswith("ACT_"):
        return "activity_background"
    if stem.startswith("SHP_TYPE"):
        return "source_geometry"
    if stem.startswith("AREA_SHP"):
        return "zoning_area"
    if stem.startswith("RATIO_INTER_INTRA"):
        return "inter_intra_ratio"
    if stem == "ATTENUATION_FORMULA":
        return "attenuation_formula"
    if stem.startswith("PLATE_SHP"):
        return "plate_shape"
    return "other"


def download(url: str, dest: Path, force: bool = False) -> None:
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "jshis-psha-input-audit/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, dest.open("wb") as out:
        shutil.copyfileobj(response, out)


def ensure_inputs(force_download: bool = False) -> tuple[Path, Path]:
    csv_zip = DOWNLOAD_DIR / CSV_ZIP
    shape_zip = DOWNLOAD_DIR / SHAPE_ZIP
    download(CSV_URL, csv_zip, force=force_download)
    download(SHAPE_URL, shape_zip, force=force_download)
    return csv_zip, shape_zip


def extract_if_needed(csv_zip: Path, shape_zip: Path, force_extract: bool = False) -> tuple[Path, Path]:
    prm_dir = EXTRACT_DIR / "prm"
    shape_dir = EXTRACT_DIR / "shape"
    if force_extract and EXTRACT_DIR.exists():
        shutil.rmtree(EXTRACT_DIR)
    if not (prm_dir / "PRM").exists():
        prm_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(csv_zip) as zf:
            zf.extractall(prm_dir)
    if not any(shape_dir.glob("*.shp")):
        shape_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(shape_zip) as zf:
            zf.extractall(shape_dir)
    return prm_dir / "PRM", shape_dir


def decode_bytes(data: bytes) -> tuple[str, str]:
    for encoding in ("utf-8-sig", "cp932", "shift_jis"):
        try:
            return data.decode(encoding), encoding
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace"), "utf-8-replace"


def csv_profile(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    text, encoding = decode_bytes(data)
    lines = text.splitlines()
    comment_lines = [line.strip() for line in lines if line.strip().startswith("#")]
    data_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith("#")]
    header = ""
    for line in comment_lines:
        clean = line.lstrip("#").strip()
        if "," in clean and not clean.upper().startswith(("UPDATED", "VER.", "DATE", "EPOCH")):
            header = clean
    if not header and data_lines:
        try:
            row = next(csv.reader([data_lines[0]]))
            header = f"{len(row)} unnamed/block fields"
        except csv.Error:
            header = "block structured"
    return {
        "encoding": encoding,
        "n_comment_lines": len(comment_lines),
        "n_data_lines": len(data_lines),
        "detected_header": header,
        "first_data_field_count": len(next(csv.reader([data_lines[0]]))) if data_lines else 0,
    }


def build_inventory(prm_dir: Path, shape_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for path in sorted(prm_dir.glob("*.csv")):
        profile = csv_profile(path)
        category_key = category_for(path.name)
        rows.append(
            {
                "file_name": path.name,
                "category": category_key,
                "category_label": CATEGORIES[category_key].label,
                "size_bytes": path.stat().st_size,
                "encoding": profile["encoding"],
                "n_comment_lines": profile["n_comment_lines"],
                "n_data_lines": profile["n_data_lines"],
                "detected_header": profile["detected_header"],
                "first_data_field_count": profile["first_data_field_count"],
            }
        )
    inventory = pd.DataFrame(rows)

    shape_rows = []
    for path in sorted(shape_dir.glob("*")):
        if path.is_file():
            shape_rows.append(
                {
                    "file_name": path.name,
                    "extension": path.suffix.lower().lstrip("."),
                    "source_layer": path.stem,
                    "size_bytes": path.stat().st_size,
                }
            )
    shape_inventory = pd.DataFrame(shape_rows)
    return inventory, shape_inventory


def build_summary(inventory: pd.DataFrame, shape_inventory: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = inventory.groupby("category", sort=False)
    for category_key, group in grouped:
        meta = CATEGORIES[category_key]
        rows.append(
            {
                "category": category_key,
                "label": meta.label,
                "n_csv_files": int(len(group)),
                "total_csv_mb": round(float(group["size_bytes"].sum()) / 1_000_000.0, 3),
                "example_file": group.iloc[0]["file_name"],
                "fields_or_block_structure": meta.fields,
                "official_role": meta.role,
                "reconstruction_status": meta.status,
                "paper_use": meta.paper_use,
            }
        )
    seen = {row["category"] for row in rows}
    for key, meta in CATEGORIES.items():
        if key not in seen and key != "other":
            rows.append(
                {
                    "category": key,
                    "label": meta.label,
                    "n_csv_files": 0,
                    "total_csv_mb": 0.0,
                    "example_file": "",
                    "fields_or_block_structure": meta.fields,
                    "official_role": meta.role,
                    "reconstruction_status": meta.status,
                    "paper_use": meta.paper_use,
                }
            )
    summary = pd.DataFrame(rows)
    shape_counts = shape_inventory["extension"].value_counts().to_dict()
    summary.attrs["shape_counts"] = shape_counts
    return summary.sort_values(["n_csv_files", "category"], ascending=[False, True]).reset_index(drop=True)


def write_download_manifest(csv_zip: Path, shape_zip: Path) -> pd.DataFrame:
    rows = [
        {
            "asset": CSV_ZIP,
            "url": CSV_URL,
            "local_cache_path": str(csv_zip.relative_to(ROOT)),
            "size_bytes": csv_zip.stat().st_size,
            "redistributed": "no",
            "use_in_this_study": "Local official parameter audit only; derived inventory tables are redistributed.",
        },
        {
            "asset": SHAPE_ZIP,
            "url": SHAPE_URL,
            "local_cache_path": str(shape_zip.relative_to(ROOT)),
            "size_bytes": shape_zip.stat().st_size,
            "redistributed": "no",
            "use_in_this_study": "Local official geometry audit only; derived inventory tables are redistributed.",
        },
        {
            "asset": "J-SHIS file format specification",
            "url": DATA_RULE_URL,
            "local_cache_path": "work/jshis_full_psha_audit/docs/A-RULES_en.pdf",
            "size_bytes": "",
            "redistributed": "no",
            "use_in_this_study": "Field definitions and file naming rules.",
        },
    ]
    return pd.DataFrame(rows)


def markdown_table(df: pd.DataFrame, columns: list[str]) -> str:
    slim = df[columns].copy()
    slim = slim.astype(str)
    out = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in slim.itertuples(index=False):
        out.append("| " + " | ".join(str(value).replace("\n", " ") for value in row) + " |")
    return "\n".join(out)


def write_report(
    summary: pd.DataFrame,
    inventory: pd.DataFrame,
    shape_inventory: pd.DataFrame,
    manifest: pd.DataFrame,
    out_path: Path,
) -> None:
    shape_counts = shape_inventory["extension"].value_counts().sort_index()
    n_shp = int(shape_counts.get("shp", 0))
    n_shape_files = int(len(shape_inventory))
    n_csv = int(len(inventory))
    total_csv_mb = inventory["size_bytes"].sum() / 1_000_000.0
    total_shape_mb = shape_inventory["size_bytes"].sum() / 1_000_000.0

    lines = [
        "# Official J-SHIS PSHA Input-Chain Audit",
        "",
        f"Audit date: 2026-06-14. Official model version: {VERSION}.",
        "",
        "## What Was Checked",
        "",
        "The audit downloads the official J-SHIS Probabilistic Seismic Hazard Maps parameter packages and records their file roles, field structures, and reconstruction status. The raw J-SHIS zip files stay in the local `work/` cache and are not redistributed in the manuscript package.",
        "",
        f"- CSV parameter package: `{CSV_ZIP}` from {CSV_URL}",
        f"- Shapefile parameter package: `{SHAPE_ZIP}` from {SHAPE_URL}",
        f"- J-SHIS download page: {DOWNLOAD_PAGE_URL}",
        f"- J-SHIS API list: {API_LIST_URL}",
        f"- File-format specification: {DATA_RULE_URL}",
        "",
        "## Main Audit Result",
        "",
        f"The official 2024 package contains {n_csv:,} CSV parameter files ({total_csv_mb:.2f} MB) and {n_shp:,} shapefile layers ({n_shape_files:,} component files, {total_shape_mb:.2f} MB). The CSV package includes source activity parameters, source geometries, occurrence-frequency grids for sources without specified faults, zoning polygons, interplate/intraplate ratios, attenuation-relation parameter codes, and Pacific/Philippine Sea plate geometry.",
        "",
        "This is stronger evidence than using only J-SHIS hazard-curve API outputs. It confirms that major official PSHA input layers are publicly accessible and versioned.",
        "",
        "The audit still does not constitute an official full J-SHIS PSHA rerun. A full rerun would require converting these J-SHIS-specific parameter files into an executable PSHA source model, validating the attenuation-relation implementation and calculation switches, and reproducing the official probability-case aggregation for all national grid cells. The manuscript should describe the present addition as an official input-chain audit and bounded reconstruction check, not as a completed official national hazard-map recomputation.",
        "",
        "## Category Summary",
        "",
        markdown_table(
            summary,
            [
                "category",
                "n_csv_files",
                "example_file",
                "fields_or_block_structure",
                "reconstruction_status",
                "paper_use",
            ],
        ),
        "",
        "## Download Manifest",
        "",
        markdown_table(manifest, ["asset", "url", "local_cache_path", "redistributed", "use_in_this_study"]),
        "",
        "## Recommended Manuscript Wording",
        "",
        "We audited the official J-SHIS 2024 PSHM parameter packages to check whether the station-residual sensitivity analysis could be connected to the official PSHA input chain. The packages provide downloadable source activity parameters, source-geometry files, occurrence-frequency grids for earthquakes without specified source faults, zoning polygons, interplate/intraplate ratios, attenuation-relation parameter codes, and Pacific/Philippine Sea plate geometry. This audit supports the provenance of the official hazard-facing check and identifies the remaining step needed for a full official rerun: translation of the J-SHIS parameter conventions into a validated executable PSHA model with the official calculation settings. The present manuscript keeps the hazard calculation at the level directly supported by public products: official response-spectrum ordinates, official hazard-curve API references, official parameter-chain provenance, and an independent OpenQuake/GEM Japan-model sensitivity run.",
        "",
        "## Files Written",
        "",
        "- `jshis_official_psha_input_audit_summary.csv`",
        "- `jshis_official_psha_input_file_inventory.csv`",
        "- `jshis_official_psha_input_shape_inventory.csv`",
        "- `jshis_official_psha_input_download_manifest.csv`",
        "- `figures/jshis_official_psha_input_chain_audit.png` and `.pdf`",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def draw_figure(summary: pd.DataFrame, fig_path_png: Path, fig_path_pdf: Path) -> None:
    plot_order = [
        "source_geometry",
        "activity_background",
        "activity_characteristic",
        "zoning_area",
        "inter_intra_ratio",
        "plate_shape",
        "eqthr",
        "attenuation_formula",
    ]
    plot_df = summary.set_index("category").reindex(plot_order).dropna(subset=["n_csv_files"]).reset_index()

    fig = plt.figure(figsize=(12.5, 7.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 0.62], wspace=0.28, hspace=0.36)
    ax_flow = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_status = fig.add_subplot(gs[1, :])

    ax_flow.axis("off")
    boxes = [
        (0.02, 0.70, 0.30, 0.20, "Official J-SHIS products\n2024 PSHM parameter zip\n2024 PSHM shapefile zip"),
        (0.36, 0.70, 0.30, 0.20, "Input layers checked\nactivity, geometry,\nfrequency, zoning,\nplate depth, GMPE codes"),
        (0.70, 0.70, 0.28, 0.20, "Derived outputs\ninventory tables\nprovenance report\nsupplementary figure"),
        (0.19, 0.30, 0.30, 0.20, "Current manuscript use\nofficial spectrum ordinates\nhazard-curve API references\nstation multiplier sensitivity"),
        (0.55, 0.30, 0.34, 0.20, "Full rerun still needs\nexecutable source conversion\nvalidated GMPE switches\nprobability-case settings"),
    ]
    for x, y, w, h, text in boxes:
        ax_flow.add_patch(plt.Rectangle((x, y), w, h, fill=False, linewidth=1.2, color="#333333"))
        ax_flow.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9.5, linespacing=1.35)
    arrows = [
        ((0.32, 0.80), (0.36, 0.80)),
        ((0.66, 0.80), (0.70, 0.80)),
        ((0.50, 0.70), (0.34, 0.50)),
        ((0.50, 0.70), (0.72, 0.50)),
    ]
    for start, end in arrows:
        ax_flow.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#333333"))
    ax_flow.text(0.0, 0.98, "A Official-parameter audit workflow", fontsize=12, fontweight="bold", va="top")

    labels = [
        "Geometry",
        "Background\nfrequency",
        "Specified-source\nactivity",
        "Zoning",
        "Inter/intra\nratio",
        "Plate shape",
        "EQTHR",
        "GMPE codes",
    ]
    values = plot_df["n_csv_files"].astype(int).to_numpy()
    y = list(range(len(values)))
    ax_bar.barh(y, values, color="#4C78A8")
    ax_bar.set_yticks(y, labels)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("CSV files")
    ax_bar.set_title("B Public J-SHIS 2024 PSHM inputs")
    ax_bar.grid(axis="x", alpha=0.25)
    for idx, value in enumerate(values):
        ax_bar.text(value + max(values) * 0.012, idx, str(int(value)), va="center", fontsize=9)

    ax_status.axis("off")
    status_text = (
        "Audit interpretation: J-SHIS publishes national PSHM parameter packages with source activity, "
        "source geometry, occurrence-frequency grids, zoning, plate-depth surfaces and attenuation-code mappings. "
        "These files support an official input-chain audit. They do not by themselves make the present paper "
        "an official full J-SHIS hazard-map recomputation, because an executable PSHA conversion and verified "
        "calculation settings are still required."
    )
    ax_status.text(0.01, 0.76, "C Reproducibility boundary", fontsize=12, fontweight="bold", va="top")
    ax_status.text(
        0.01,
        0.46,
        status_text,
        fontsize=10.5,
        va="center",
        ha="left",
        wrap=True,
        bbox=dict(boxstyle="square,pad=0.55", facecolor="#F5F5F5", edgecolor="#333333", linewidth=1.0),
    )

    fig.savefig(fig_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(fig_path_pdf, bbox_inches="tight")
    plt.close(fig)


def run(force_download: bool = False, force_extract: bool = False) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    csv_zip, shape_zip = ensure_inputs(force_download=force_download)
    prm_dir, shape_dir = extract_if_needed(csv_zip, shape_zip, force_extract=force_extract)
    inventory, shape_inventory = build_inventory(prm_dir, shape_dir)
    summary = build_summary(inventory, shape_inventory)
    manifest = write_download_manifest(csv_zip, shape_zip)

    summary.to_csv(OUTPUT_DIR / "jshis_official_psha_input_audit_summary.csv", index=False)
    inventory.to_csv(OUTPUT_DIR / "jshis_official_psha_input_file_inventory.csv", index=False)
    shape_inventory.to_csv(OUTPUT_DIR / "jshis_official_psha_input_shape_inventory.csv", index=False)
    manifest.to_csv(OUTPUT_DIR / "jshis_official_psha_input_download_manifest.csv", index=False)
    write_report(
        summary,
        inventory,
        shape_inventory,
        manifest,
        OUTPUT_DIR / "jshis_official_psha_input_chain_audit.md",
    )
    draw_figure(
        summary,
        FIGURE_DIR / "jshis_official_psha_input_chain_audit.png",
        FIGURE_DIR / "jshis_official_psha_input_chain_audit.pdf",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force-download", action="store_true", help="Download official J-SHIS packages even if cached.")
    parser.add_argument("--force-extract", action="store_true", help="Recreate extracted working copies.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(force_download=args.force_download, force_extract=args.force_extract)


if __name__ == "__main__":
    main()
