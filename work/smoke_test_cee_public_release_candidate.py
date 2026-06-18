#!/usr/bin/env python3
"""Smoke-test the CEE public release candidate v0.8 package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import py_compile
import shutil
import subprocess
import zipfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ZIP = ROOT / "outputs" / "cee_public_release_candidate_v0_8.zip"
DEFAULT_EXTRACT_DIR = ROOT / "work" / "cee_public_release_candidate_v0_8_smoke"
OUT_CSV = ROOT / "outputs" / "cee_public_release_candidate_v0_8_smoke_test.csv"
OUT_MD = ROOT / "outputs" / "cee_public_release_candidate_v0_8_smoke_test.md"

ARTICLE_DIR = Path("outputs/cee_submission_latex_v0_8_english_article")
FORBIDDEN_FILE_PATTERNS = [
    ".mseed",
    ".sac",
    ".BSON",
    "__pycache__",
    ".pyc",
    "flatfile_sub1-v2024.zip",
    "flatfile-v2024.zip",
    "flatfile_sub2-v2024.zip",
    "smrec_schema.tsv",
    "site_schema.tsv",
    "source_schema.tsv",
    "gsi_dem_png_z",
    "external_data/",
    "extracted/",
]
FORBIDDEN_TEXT_TOKENS = ["/Users/", "/private/var/", "/var/folders/"]


@dataclass
class Check:
    name: str
    status: str
    detail: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add(checks: list[Check], name: str, ok: bool, detail: str) -> None:
    checks.append(Check(name, "PASS" if ok else "FAIL", detail.replace("\n", " ")))


def extract_zip(zip_path: Path, extract_dir: Path) -> None:
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(extract_dir)


def pdf_pages(path: Path) -> int | None:
    try:
        output = subprocess.check_output(["pdfinfo", str(path)], text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    for line in output.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    return None


def iter_text_files(base: Path) -> list[Path]:
    suffixes = {".md", ".py", ".txt", ".tex", ".csv", ".json", ".yml", ".yaml"}
    return [path for path in base.rglob("*") if path.is_file() and path.suffix.lower() in suffixes]


def scan_forbidden_files(release_dir: Path, checks: list[Check]) -> None:
    files = [path for path in release_dir.rglob("*") if path.is_file()]
    filename_hits: list[str] = []
    text_hits: list[str] = []
    nested_zip_hits: list[str] = []
    nested_text_hits: list[str] = []

    for path in files:
        rel = str(path.relative_to(release_dir))
        lower = rel.lower()
        if any(pattern.lower() in lower for pattern in FORBIDDEN_FILE_PATTERNS):
            filename_hits.append(rel)

    for path in iter_text_files(release_dir):
        text = path.read_text(encoding="utf-8", errors="ignore")
        if any(token in text for token in FORBIDDEN_TEXT_TOKENS):
            text_hits.append(str(path.relative_to(release_dir)))

    for path in files:
        if path.suffix.lower() != ".zip":
            continue
        with zipfile.ZipFile(path) as archive:
            bad_member = archive.testzip()
            if bad_member:
                nested_zip_hits.append(f"{path.name}: bad CRC {bad_member}")
            for info in archive.infolist():
                name = info.filename
                lower_name = name.lower()
                if any(pattern.lower() in lower_name for pattern in FORBIDDEN_FILE_PATTERNS):
                    nested_zip_hits.append(f"{path.name}: {name}")
                if Path(name).suffix.lower() in {".md", ".py", ".txt", ".tex", ".csv", ".json", ".yml", ".yaml"}:
                    data = archive.read(info)
                    text = data.decode("utf-8", errors="ignore")
                    if any(token in text for token in FORBIDDEN_TEXT_TOKENS):
                        nested_text_hits.append(f"{path.name}: {name}")

    add(checks, "forbidden raw/cache filenames", not filename_hits, "none" if not filename_hits else ", ".join(filename_hits[:10]))
    add(checks, "local absolute paths", not text_hits, "none" if not text_hits else ", ".join(text_hits[:10]))
    add(checks, "nested zip raw/cache scan", not nested_zip_hits, "none" if not nested_zip_hits else ", ".join(nested_zip_hits[:10]))
    add(checks, "nested zip local paths", not nested_text_hits, "none" if not nested_text_hits else ", ".join(nested_text_hits[:10]))


def expected_file_checks(release_dir: Path, checks: list[Check]) -> None:
    expected = [
        "README.md",
        "outputs/cee_submission_latex_v0_8_english_article.zip",
        "outputs/cee_submission_package_sha256.txt",
        "outputs/cee_v08_public_archive_readiness.md",
        "outputs/cee_v08_public_archive_file_review.csv",
        str(ARTICLE_DIR / "README.md"),
        str(ARTICLE_DIR / "main.tex"),
        str(ARTICLE_DIR / "main.pdf"),
        str(ARTICLE_DIR / "supplementary_information.tex"),
        str(ARTICLE_DIR / "supplementary_information.pdf"),
        str(ARTICLE_DIR / "cover_letter_cee.md"),
        str(ARTICLE_DIR / "build_figure8_en.py"),
        str(ARTICLE_DIR / "build_figure9_en.py"),
        str(ARTICLE_DIR / "figures/figure8.pdf"),
        str(ARTICLE_DIR / "figures/figure8.png"),
        str(ARTICLE_DIR / "figures/figure9_engineering_cases_uncertainty.pdf"),
        str(ARTICLE_DIR / "figures/figure9_engineering_cases_uncertainty.png"),
        str(ARTICLE_DIR / "supplement/cee_claim_evidence_traceability.csv"),
        str(ARTICLE_DIR / "supplement/kiknet_surface_downhole_ratios.csv"),
        str(ARTICLE_DIR / "supplement/kiknet_surface_downhole_audit.md"),
        str(ARTICLE_DIR / "supplement/kiknet_surface_downhole_audit.py"),
        str(ARTICLE_DIR / "supplement/kiknet_multievent_transfer_functions.csv"),
        str(ARTICLE_DIR / "supplement/kiknet_multievent_transfer_function_bins.csv"),
        str(ARTICLE_DIR / "supplement/kiknet_multievent_transfer_function_event_summary.csv"),
        str(ARTICLE_DIR / "supplement/kiknet_multievent_transfer_function_station_summary.csv"),
        str(ARTICLE_DIR / "supplement/kiknet_multievent_transfer_function_audit.md"),
        str(ARTICLE_DIR / "supplement/kiknet_multievent_transfer_function_audit.py"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_nonergodic_psha_summary.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_station_corrections.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_station_surface_validation.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_station_surface_validation_summary.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_continuous_correction_surface_mesh.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_official_response_station_values.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_official_response_mesh_values.csv"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_nonergodic_psha_audit.md"),
        str(ARTICLE_DIR / "supplement/jshis_multiperiod_nonergodic_psha_audit.py"),
        str(ARTICLE_DIR / "supplement/jshis_official_psha_input_chain_audit.md"),
        str(ARTICLE_DIR / "supplement/jshis_official_psha_input_audit_summary.csv"),
        str(ARTICLE_DIR / "supplement/jshis_official_psha_input_file_inventory.csv"),
        str(ARTICLE_DIR / "supplement/jshis_official_psha_input_shape_inventory.csv"),
        str(ARTICLE_DIR / "supplement/jshis_official_psha_input_download_manifest.csv"),
        str(ARTICLE_DIR / "supplement/jshis_official_psha_input_audit.py"),
        str(ARTICLE_DIR / "supplement/jshis_representative_city_spectrum_cases.csv"),
        str(ARTICLE_DIR / "supplement/jshis_representative_uncertainty_summary.csv"),
    ]
    missing = [item for item in expected if not (release_dir / item).exists()]
    add(checks, "expected files", not missing, f"{len(expected)} expected files present" if not missing else ", ".join(missing))


def article_zip_checks(release_dir: Path, checks: list[Check]) -> None:
    zip_path = release_dir / "outputs" / "cee_submission_latex_v0_8_english_article.zip"
    sha_path = release_dir / "outputs" / "cee_submission_package_sha256.txt"
    if not zip_path.exists() or not sha_path.exists():
        add(checks, "article zip checksum", False, "zip or checksum file missing")
        return
    recorded = sha_path.read_text(encoding="utf-8").split()[0]
    actual = sha256(zip_path)
    add(checks, "article zip checksum", recorded == actual, actual)
    with zipfile.ZipFile(zip_path) as archive:
        names = set(archive.namelist())
        required = {
            "main.pdf",
            "supplementary_information.pdf",
            "figures/figure8.pdf",
            "figures/figure9_engineering_cases_uncertainty.pdf",
            "supplement/kiknet_surface_downhole_ratios.csv",
            "supplement/kiknet_surface_downhole_audit.py",
            "supplement/kiknet_multievent_transfer_functions.csv",
            "supplement/kiknet_multievent_transfer_function_audit.py",
            "supplement/jshis_multiperiod_nonergodic_psha_summary.csv",
            "supplement/jshis_multiperiod_nonergodic_psha_audit.py",
            "supplement/jshis_official_psha_input_chain_audit.md",
            "supplement/jshis_official_psha_input_audit_summary.csv",
            "supplement/jshis_official_psha_input_audit.py",
            "supplement/jshis_representative_city_spectrum_cases.csv",
            "supplement/jshis_representative_uncertainty_summary.csv",
        }
        missing = sorted(required - names)
        add(checks, "article zip members", not missing, "required members present" if not missing else ", ".join(missing))


def pdf_checks(release_dir: Path, checks: list[Check]) -> None:
    main_pdf = release_dir / ARTICLE_DIR / "main.pdf"
    supp_pdf = release_dir / ARTICLE_DIR / "supplementary_information.pdf"
    main_pages = pdf_pages(main_pdf)
    supp_pages = pdf_pages(supp_pdf)
    add(checks, "PDF page counts", main_pages == 17 and supp_pages == 14, f"main={main_pages}; supplementary={supp_pages}")


def table_checks(release_dir: Path, checks: list[Check]) -> None:
    supplement = release_dir / ARTICLE_DIR / "supplement"
    failures: list[str] = []

    trace = pd.read_csv(supplement / "cee_claim_evidence_traceability.csv")
    if trace.shape != (17, 7):
        failures.append(f"traceability shape {trace.shape}")
    if {"C13", "C14", "C15", "C16", "C17"} - set(trace["claim_id"].astype(str)):
        failures.append("traceability missing C13/C14/C15/C16/C17")

    kiknet = pd.read_csv(supplement / "kiknet_surface_downhole_ratios.csv")
    if kiknet.shape[0] != 288:
        failures.append(f"kiknet rows {kiknet.shape[0]}")
    for col in ["spectral_ratio_0.2_0.5_hz", "spectral_ratio_0.5_1_hz", "horizontal_pga_ratio"]:
        if col not in kiknet.columns:
            failures.append(f"kiknet missing {col}")
    if "spectral_ratio_0.2_0.5_hz" in kiknet.columns:
        med = float(kiknet["spectral_ratio_0.2_0.5_hz"].median())
        if abs(med - 2.382) > 0.01:
            failures.append(f"kiknet 0.2-0.5 Hz median {med:.3f}")

    kiknet_multi = pd.read_csv(supplement / "kiknet_multievent_transfer_functions.csv")
    if kiknet_multi.shape[0] != 460:
        failures.append(f"kiknet multi-event rows {kiknet_multi.shape[0]}")
    if kiknet_multi["event_id"].nunique() != 2:
        failures.append(f"kiknet multi-event count {kiknet_multi['event_id'].nunique()}")
    med_multi = float(kiknet_multi["horizontal_pga_ratio"].median())
    if abs(med_multi - 1.458) > 0.01:
        failures.append(f"kiknet multi-event horizontal PGA median {med_multi:.3f}")

    summary = pd.read_csv(supplement / "jshis_official_hazard_response_summary.csv")
    if summary.shape != (6, 11):
        failures.append(f"response summary shape {summary.shape}")

    multiperiod = pd.read_csv(supplement / "jshis_multiperiod_nonergodic_psha_summary.csv")
    grid_50y10 = multiperiod[
        multiperiod["scope"].eq("sampled_grid")
        & multiperiod["probability_level"].eq("50y_10pct")
        & multiperiod["period_s"].notna()
    ]
    if grid_50y10.shape[0] != 8:
        failures.append(f"multiperiod sampled-grid 50y10 rows {grid_50y10.shape[0]}")
    for period, expected in [(0.3, 1.045), (1.0, 0.834), (3.0, 0.658), (5.0, 0.658)]:
        sub = grid_50y10[grid_50y10["period_s"].round(3).eq(period)]
        if sub.empty:
            failures.append(f"multiperiod missing {period}s")
        elif abs(float(sub.iloc[0]["multiplier_q50"]) - expected) > 0.01:
            failures.append(f"multiperiod {period}s multiplier {float(sub.iloc[0]['multiplier_q50']):.3f}")

    input_summary = pd.read_csv(supplement / "jshis_official_psha_input_audit_summary.csv")
    if input_summary.shape[0] != 8:
        failures.append(f"official input-chain categories {input_summary.shape[0]}")
    if int(input_summary["n_csv_files"].sum()) != 421:
        failures.append(f"official input-chain CSV files {int(input_summary['n_csv_files'].sum())}")
    shape_inventory = pd.read_csv(supplement / "jshis_official_psha_input_shape_inventory.csv")
    if shape_inventory.shape[0] != 800:
        failures.append(f"official input-chain shapefile components {shape_inventory.shape[0]}")

    cases = pd.read_csv(supplement / "jshis_representative_city_spectrum_cases.csv")
    if cases.shape[0] != 40 or cases["case_name"].nunique() != 5:
        failures.append(f"representative cases shape {cases.shape}; cases {cases['case_name'].nunique()}")
    tokyo_sa3 = cases[
        cases["case_name"].eq("Tokyo") & cases["period_s"].round(3).eq(3.0)
    ]
    if tokyo_sa3.empty or abs(float(tokyo_sa3.iloc[0]["delta_pct_at_ordinate"]) + 34.64) > 0.1:
        failures.append("representative Tokyo SA3 delta mismatch")

    uncertainty = pd.read_csv(supplement / "jshis_representative_uncertainty_summary.csv")
    med_delta = uncertainty[
        uncertainty["quantity"].eq("delta_pct") & uncertainty["quantile"].round(2).eq(0.50)
    ]
    if med_delta.empty or abs(float(med_delta.iloc[0]["observed"]) + 38.09) > 0.1:
        failures.append("representative uncertainty median delta mismatch")

    add(checks, "table schemas and key values", not failures, "traceability, KiK-net, and response-spectrum tables verified" if not failures else "; ".join(failures))


def script_checks(release_dir: Path, checks: list[Check]) -> None:
    scripts = [
        release_dir / ARTICLE_DIR / "build_figure8_en.py",
        release_dir / ARTICLE_DIR / "build_figure9_en.py",
        release_dir / ARTICLE_DIR / "supplement/kiknet_surface_downhole_audit.py",
        release_dir / ARTICLE_DIR / "supplement/kiknet_multievent_transfer_function_audit.py",
        release_dir / ARTICLE_DIR / "supplement/jshis_multiperiod_nonergodic_psha_audit.py",
        release_dir / ARTICLE_DIR / "supplement/jshis_official_psha_input_audit.py",
    ]
    failures = []
    for script in scripts:
        if not script.exists():
            failures.append(f"missing {script.relative_to(release_dir)}")
            continue
        try:
            py_compile.compile(str(script), doraise=True)
        except py_compile.PyCompileError as exc:
            failures.append(f"{script.name}: {exc.msg}")
    add(checks, "script syntax", not failures, f"{len(scripts)} scripts compile" if not failures else "; ".join(failures))


def write_reports(checks: list[Check], zip_path: Path, release_dir: Path) -> None:
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"])
        writer.writeheader()
        for check in checks:
            writer.writerow({"check": check.name, "status": check.status, "detail": check.detail})

    failures = [check for check in checks if check.status != "PASS"]
    lines = [
        "# CEE public release candidate v0.8 smoke test",
        "",
        f"Generated: {date.today().isoformat()}.",
        "",
        f"- Zip: `{zip_path.relative_to(ROOT)}`.",
        f"- Extracted test directory: `{release_dir.relative_to(ROOT)}`.",
        f"- Status: {'passed' if not failures else 'failed'}.",
        "",
        "| Check | Status | Detail |",
        "| --- | ---: | --- |",
    ]
    for check in checks:
        lines.append(f"| {check.name} | {check.status} | {check.detail} |")
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "This smoke test checks package structure, key manuscript PDFs, critical derived tables, nested zip integrity, and absence of raw waveform/cache files or local absolute paths. It does not rerun the scientific workflows.",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def run(zip_path: Path, extract_dir: Path) -> tuple[list[Check], Path]:
    checks: list[Check] = []
    add(checks, "zip exists", zip_path.exists(), zip_path.name)
    if not zip_path.exists():
        return checks, extract_dir
    add(checks, "zip is readable", zipfile.is_zipfile(zip_path), "zipfile header valid")
    if not zipfile.is_zipfile(zip_path):
        return checks, extract_dir
    with zipfile.ZipFile(zip_path) as archive:
        bad = archive.testzip()
    add(checks, "zip member CRC", bad is None, "all members OK" if bad is None else f"bad member: {bad}")
    extract_zip(zip_path, extract_dir)
    add(checks, "release directory extracted", extract_dir.exists(), extract_dir.name)
    expected_file_checks(extract_dir, checks)
    scan_forbidden_files(extract_dir, checks)
    article_zip_checks(extract_dir, checks)
    pdf_checks(extract_dir, checks)
    table_checks(extract_dir, checks)
    script_checks(extract_dir, checks)
    return checks, extract_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Smoke-test the CEE v0.8 public release candidate.")
    parser.add_argument("--zip-path", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--extract-dir", type=Path, default=DEFAULT_EXTRACT_DIR)
    args = parser.parse_args()
    checks, release_dir = run(args.zip_path.resolve(), args.extract_dir.resolve())
    write_reports(checks, args.zip_path.resolve(), release_dir)
    print(OUT_MD.relative_to(ROOT))
    if any(check.status != "PASS" for check in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
