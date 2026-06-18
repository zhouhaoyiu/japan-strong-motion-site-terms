#!/usr/bin/env python3
"""Build and audit a curated public-release candidate from Git-visible files."""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
import zipfile
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "outputs"
OUT_CSV = OUT_DIR / "cee_public_release_candidate_v0_8_audit.csv"
OUT_MD = OUT_DIR / "cee_public_release_candidate_v0_8_audit.md"
OUT_ZIP = OUT_DIR / "cee_public_release_candidate_v0_8.zip"
OUT_SHA = OUT_DIR / "cee_public_release_candidate_v0_8_sha256.txt"

GENERATED = {
    "outputs/cee_public_release_candidate_v0_8_audit.csv",
    "outputs/cee_public_release_candidate_v0_8_audit.md",
    "outputs/cee_public_release_candidate_v0_8_smoke_test.csv",
    "outputs/cee_public_release_candidate_v0_8_smoke_test.md",
    "outputs/cee_public_release_candidate_v0_8.zip",
    "outputs/cee_public_release_candidate_v0_8_sha256.txt",
}

RAW_PATTERNS = [
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


def git_visible_files() -> list[str]:
    out = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        text=True,
    )
    files = []
    for line in out.splitlines():
        if not line or line in GENERATED:
            continue
        path = ROOT / line
        if path.is_file():
            files.append(line)
    return sorted(files)


def decode_text(data: bytes) -> str | None:
    if b"\x00" in data[:4096]:
        return None
    for encoding in ("utf-8", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return None


def text_flags(path: Path) -> tuple[str, str]:
    text = decode_text(path.read_bytes())
    if text is None:
        return "", ""
    abs_hits = sorted({pat.pattern for pat in ABSOLUTE_PATH_PATTERNS if pat.search(text)})
    sensitive_hits = sorted({pat.pattern for pat in SENSITIVE_TEXT_PATTERNS if pat.search(text)})
    return ";".join(abs_hits), ";".join(sensitive_hits)


def raw_flag(rel: str) -> str:
    return "yes" if any(pattern.search(rel) for pattern in RAW_PATTERNS) else "no"


def zip_flags(path: Path) -> tuple[int, str]:
    if path.suffix.lower() != ".zip":
        return 0, ""
    bad_entries: list[str] = []
    with zipfile.ZipFile(path) as zf:
        names = [info.filename for info in zf.infolist()]
    for name in names:
        if name.startswith("outputs/") or name.startswith("/") or raw_flag(name) == "yes":
            bad_entries.append(name)
        if any(token in name for token in ("cee_80_acceptance", "reviewer_risk", "v08_80", "chinese_review")):
            bad_entries.append(name)
    return len(names), ";".join(sorted(set(bad_entries))[:20])


def category(rel: str) -> str:
    lower = rel.lower()
    if lower.endswith(".py"):
        return "script"
    if lower.endswith((".tex", ".md", ".json", ".yml", ".yaml", ".txt")):
        return "source_or_documentation"
    if lower.endswith(".csv"):
        return "derived_table"
    if lower.endswith((".pdf", ".png")):
        return "rendered_output"
    if lower.endswith(".zip"):
        return "package_zip"
    return "other"


def build_rows(files: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rel in files:
        path = ROOT / rel
        abs_hits, sensitive_hits = text_flags(path)
        zip_entry_count, zip_bad_entries = zip_flags(path)
        rows.append(
            {
                "path": rel,
                "category": category(rel),
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256_file(path),
                "raw_or_cache_flag": raw_flag(rel),
                "absolute_path_patterns": abs_hits,
                "sensitive_text_patterns": sensitive_hits,
                "zip_entry_count": str(zip_entry_count),
                "zip_bad_entries": zip_bad_entries,
            }
        )
    return rows


def write_csv(rows: list[dict[str, str]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "path",
        "category",
        "size_bytes",
        "sha256",
        "raw_or_cache_flag",
        "absolute_path_patterns",
        "sensitive_text_patterns",
        "zip_entry_count",
        "zip_bad_entries",
    ]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_zip(files: list[str]) -> str:
    if OUT_ZIP.exists():
        OUT_ZIP.unlink()
    with zipfile.ZipFile(OUT_ZIP, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for rel in files:
            zf.write(ROOT / rel, rel)
    digest = sha256_file(OUT_ZIP)
    OUT_SHA.write_text(f"{digest}  {OUT_ZIP.relative_to(ROOT)}\n", encoding="utf-8")
    return digest


def write_md(rows: list[dict[str, str]], digest: str) -> None:
    raw_hits = [row for row in rows if row["raw_or_cache_flag"] == "yes"]
    abs_hits = [row for row in rows if row["absolute_path_patterns"]]
    sensitive_hits = [row for row in rows if row["sensitive_text_patterns"]]
    zip_bad = [row for row in rows if row["zip_bad_entries"]]
    total_size = sum(int(row["size_bytes"]) for row in rows)
    text = [
        "# CEE public release candidate v0.8 audit",
        "",
        f"Generated: {date.today().isoformat()}.",
        "",
        "## Summary",
        "",
        f"- Candidate files: {len(rows)}.",
        f"- Total file bytes before zip: {total_size}.",
        f"- Local candidate zip: `{OUT_ZIP.relative_to(ROOT)}`.",
        f"- Local candidate zip SHA-256: `{digest}`.",
        f"- Raw/cache path flags: {len(raw_hits)}.",
        f"- Absolute local path flags: {len(abs_hits)}.",
        f"- Private-planning text flags: {len(sensitive_hits)}.",
        f"- Zip member-name flags: {len(zip_bad)}.",
        "",
        "## Interpretation",
        "",
        "This audit is built from Git-visible files after applying `.gitignore`. It checks the candidate public repository surface, not the untracked local working files. A clean result supports GitHub/Zenodo preparation, but final author/legal approval is still required for derived tables.",
        "",
    ]
    if raw_hits or abs_hits or sensitive_hits or zip_bad:
        text.extend(["## Open flags", ""])
        for label, hits, field in [
            ("raw/cache", raw_hits, "raw_or_cache_flag"),
            ("absolute path", abs_hits, "absolute_path_patterns"),
            ("sensitive text", sensitive_hits, "sensitive_text_patterns"),
            ("zip member", zip_bad, "zip_bad_entries"),
        ]:
            if hits:
                text.append(f"### {label}")
                for row in hits[:30]:
                    text.append(f"- `{row['path']}`: {row[field]}")
                text.append("")
    OUT_MD.write_text("\n".join(text), encoding="utf-8")


def main() -> None:
    files = git_visible_files()
    rows = build_rows(files)
    digest = write_zip(files)
    write_csv(rows)
    write_md(rows, digest)
    print(f"wrote {OUT_CSV.relative_to(ROOT)}")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    print(f"wrote {OUT_ZIP.relative_to(ROOT)}")
    print(f"sha256 {digest}")


if __name__ == "__main__":
    main()
