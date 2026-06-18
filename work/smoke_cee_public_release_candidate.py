#!/usr/bin/env python3
"""Smoke-test the curated CEE public-release candidate zip."""

from __future__ import annotations

import csv
import hashlib
import py_compile
import re
import shutil
import subprocess
import zipfile
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_ZIP = ROOT / "outputs" / "cee_public_release_candidate_v0_8.zip"
CANDIDATE_SHA = ROOT / "outputs" / "cee_public_release_candidate_v0_8_sha256.txt"
AUDIT_MD = ROOT / "outputs" / "cee_public_release_candidate_v0_8_audit.md"
OUT_MD = ROOT / "outputs" / "cee_public_release_candidate_v0_8_smoke_test.md"
OUT_CSV = ROOT / "outputs" / "cee_public_release_candidate_v0_8_smoke_test.csv"
EXTRACT_DIR = Path("/private/tmp/cee_public_release_candidate_v0_8_smoke")

NESTED_ZIP = "outputs/cee_submission_latex_v0_8_english_article.zip"
NESTED_SHA = "outputs/cee_submission_package_sha256.txt"
MAIN_PDF = "outputs/cee_submission_latex_v0_8_english_article/main.pdf"
SUPP_PDF = "outputs/cee_submission_latex_v0_8_english_article/supplementary_information.pdf"

REQUIRED = [
    "README.md",
    ".gitignore",
    "outputs/cee_submission_latex_v0_8_english_article/main.tex",
    "outputs/cee_submission_latex_v0_8_english_article/supplementary_information.tex",
    "outputs/cee_submission_latex_v0_8_english_article/main.pdf",
    "outputs/cee_submission_latex_v0_8_english_article/supplementary_information.pdf",
    "outputs/cee_submission_latex_v0_8_english_article/cover_letter_cee.md",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/kiknet_multievent_transfer_functions.csv",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_multiperiod_nonergodic_psha_summary.csv",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_official_psha_input_chain_audit.md",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_official_psha_input_audit_summary.csv",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_official_psha_input_audit.py",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_representative_city_spectrum_cases.csv",
    "outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_representative_uncertainty_summary.csv",
    "outputs/cee_submission_latex_v0_8_english_article/figures/figure9_engineering_cases_uncertainty.pdf",
    "outputs/cee_submission_latex_v0_8_english_article/build_figure9_en.py",
    "outputs/cee_submission_latex_v0_8_english_article.zip",
    "outputs/cee_submission_package_sha256.txt",
]

RAW_CACHE_PATTERNS = [
    re.compile(pattern, re.IGNORECASE)
    for pattern in [
        r"flatfile.*\.zip$",
        r"smrec_schema\.tsv$",
        r"source_schema\.tsv$",
        r"site_schema\.tsv$",
        r"\.mseed$",
        r"\.sac$",
        r"waveform",
        r"seisbench",
        r"gsi_dem_png_z",
        r"external_data/",
        r"extracted/",
        r"__pycache__",
    ]
]

ABSOLUTE_PATH_PATTERNS = [
    re.compile(r"/Users/yojironoda/"),
    re.compile(r"/private/var/"),
    re.compile(r"/var/folders/"),
]

SENSITIVE_TEXT_PATTERNS = [
    re.compile(pattern, re.IGNORECASE)
    for pattern in [
        r"advisor approval",
        r"advisor accepted",
        r"advisor/domain",
        r"near-80",
        r"80% acceptance",
        r"acceptance[-\s]+probability",
        r"true acceptance probability",
        r"CEE true acceptance",
        r"reviewer strategy",
        r"advisor feedback",
        r"Submission probability effect",
        r"probability effect",
        r"guaranteed acceptance",
        r"导师",
        r"录用",
        r"保证录用",
        r"投稿概率",
    ]
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_expected_sha() -> str:
    first = CANDIDATE_SHA.read_text(encoding="utf-8").split()[0]
    if not re.fullmatch(r"[0-9a-f]{64}", first):
        raise ValueError(f"bad candidate sha file: {CANDIDATE_SHA}")
    return first


def add_gate(rows: list[dict[str, str]], gate: str, status: str, detail: str) -> None:
    rows.append({"gate": gate, "status": status, "detail": detail})


def pdf_pages(path: Path) -> int:
    out = subprocess.check_output(["pdfinfo", str(path)], text=True)
    for line in out.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"could not read page count from {path}")


def decode_text(data: bytes) -> str | None:
    if b"\x00" in data[:4096]:
        return None
    for encoding in ("utf-8", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return None


def text_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.stat().st_size > 3_000_000:
            continue
        if path.suffix.lower() in {".pdf", ".png", ".zip"}:
            continue
        files.append(path)
    return files


def scan_text(root: Path, patterns: list[re.Pattern[str]]) -> list[str]:
    hits: list[str] = []
    for path in text_files(root):
        text = decode_text(path.read_bytes())
        if text is None:
            continue
        rel = path.relative_to(root).as_posix()
        for pattern in patterns:
            if pattern.search(text):
                hits.append(f"{rel}:{pattern.pattern}")
    return sorted(set(hits))


def raw_cache_name_hits(names: list[str]) -> list[str]:
    return sorted({name for name in names if any(pattern.search(name) for pattern in RAW_CACHE_PATTERNS)})


def nested_zip_bad_entries(path: Path) -> tuple[int, list[str]]:
    bad: list[str] = []
    with zipfile.ZipFile(path) as zf:
        names = [info.filename for info in zf.infolist()]
    for name in names:
        lower = name.lower()
        if name.startswith("/") or name.startswith("outputs/") or "__pycache__" in lower:
            bad.append(name)
        if lower.endswith((".aux", ".log", ".out", ".fls", ".fdb_latexmk", ".synctex.gz")):
            bad.append(name)
        if any(token in lower for token in ["acceptance", "reviewer_risk", "advisor", "chinese_review"]):
            bad.append(name)
    return len(names), sorted(set(bad))


def compile_python(root: Path) -> tuple[int, list[str]]:
    failures: list[str] = []
    count = 0
    for path in root.rglob("*.py"):
        count += 1
        try:
            py_compile.compile(str(path), doraise=True)
        except py_compile.PyCompileError as exc:
            failures.append(f"{path.relative_to(root).as_posix()}: {exc.msg}")
    return count, failures


def main() -> None:
    rows: list[dict[str, str]] = []

    expected = read_expected_sha()
    actual = sha256_file(CANDIDATE_ZIP)
    add_gate(rows, "candidate_zip_checksum", "pass" if actual == expected else "fail", actual)

    if EXTRACT_DIR.exists():
        shutil.rmtree(EXTRACT_DIR)
    EXTRACT_DIR.mkdir(parents=True)
    with zipfile.ZipFile(CANDIDATE_ZIP) as zf:
        names = [info.filename for info in zf.infolist()]
        zf.extractall(EXTRACT_DIR)

    audit_text = AUDIT_MD.read_text(encoding="utf-8")
    match = re.search(r"Candidate files:\s+(\d+)", audit_text)
    audit_count = int(match.group(1)) if match else -1
    add_gate(rows, "candidate_file_count", "pass" if len(names) == audit_count else "fail", f"zip={len(names)} audit={audit_count}")

    missing = [name for name in REQUIRED if name not in names]
    add_gate(rows, "required_files", "pass" if not missing else "fail", "all present" if not missing else ",".join(missing))

    nested_sha_expected = (EXTRACT_DIR / NESTED_SHA).read_text(encoding="utf-8").split()[0]
    nested_sha_actual = sha256_file(EXTRACT_DIR / NESTED_ZIP)
    add_gate(rows, "nested_article_zip_checksum", "pass" if nested_sha_actual == nested_sha_expected else "fail", f"{NESTED_ZIP}: {nested_sha_actual}")

    nested_count, nested_bad = nested_zip_bad_entries(EXTRACT_DIR / NESTED_ZIP)
    add_gate(rows, "nested_article_zip_entries", "pass" if not nested_bad else "fail", f"entries={nested_count} bad={nested_bad[:20]}")

    main_pages = pdf_pages(EXTRACT_DIR / MAIN_PDF)
    supp_pages = pdf_pages(EXTRACT_DIR / SUPP_PDF)
    add_gate(rows, "pdf_pages", "pass" if (main_pages == 18 and supp_pages == 14) else "fail", f"main={main_pages}; supplement={supp_pages}")

    raw_hits = raw_cache_name_hits(names)
    add_gate(rows, "raw_cache_path_scan", "pass" if not raw_hits else "fail", f"hits={len(raw_hits)}")

    abs_hits = scan_text(EXTRACT_DIR, ABSOLUTE_PATH_PATTERNS)
    add_gate(rows, "absolute_path_scan", "pass" if not abs_hits else "fail", f"hits={len(abs_hits)}")

    sensitive_hits = scan_text(EXTRACT_DIR, SENSITIVE_TEXT_PATTERNS)
    add_gate(rows, "sensitive_text_scan", "pass" if not sensitive_hits else "fail", f"hits={len(sensitive_hits)}")

    py_count, py_failures = compile_python(EXTRACT_DIR)
    add_gate(rows, "python_syntax_compile", "pass" if not py_failures else "fail", f"files={py_count} failures={len(py_failures)}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gate", "status", "detail"])
        writer.writeheader()
        writer.writerows(rows)

    pass_count = sum(1 for row in rows if row["status"] == "pass")
    text = [
        "# CEE public release candidate v0.8 smoke test",
        "",
        f"Generated: {date.today().isoformat()}.",
        "",
        "## Summary",
        "",
        f"- Gates checked: {len(rows)}.",
        f"- Passing gates: {pass_count}.",
        f"- Open or review gates: {len(rows) - pass_count}.",
        f"- Tested extracted directory: `{EXTRACT_DIR}`.",
        f"- Candidate zip SHA-256: `{actual}`.",
        "",
        "## Gate table",
        "",
    ]
    for row in rows:
        text.append(f"- `{row['gate']}`: {row['status']} ({row['detail']})")
    if sensitive_hits or abs_hits:
        text.extend(["", "## Scan hits", ""])
        for hit in (abs_hits + sensitive_hits)[:80]:
            text.append(f"- `{hit}`")
    OUT_MD.write_text("\n".join(text) + "\n", encoding="utf-8")

    print(f"wrote {OUT_CSV.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    print(f"passes {pass_count}/{len(rows)}")
    if pass_count != len(rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
