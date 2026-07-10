#!/usr/bin/env python3
"""Build the compiled CEE upload files and a compact peer-review archive."""

from __future__ import annotations

import csv
import hashlib
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
OUTPUT = ROOT / "outputs" / "cee_initial_submission_2026-07-11"
PREFIX = "Zhou_Ma_Repeatable_station_terms_CEE"
MAX_DERIVED_BYTES = 2_500_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reviewer_files() -> tuple[list[Path], list[Path]]:
    included = [
        ROOT / ".github" / "workflows" / "ci.yml",
        ROOT / "LICENSE",
        ROOT / "README.md",
        ROOT / "environment.yml",
        ROOT / "public_inputs_manifest.tsv",
        ARTICLE / "README.md",
        ARTICLE / "main.tex",
        ARTICLE / "supplementary_information.tex",
        ARTICLE / "references_shared.tex",
        ARTICLE / "reference_traceability.md",
        ARTICLE / "build_event_adjusted_supplement_figures.py",
    ]
    included.extend(sorted((ROOT / "work").glob("*.py")))
    included.extend(sorted((ARTICLE / "supplement").glob("*.md")))
    included.extend(sorted((ARTICLE / "supplement").glob("*.py")))

    large_derived: list[Path] = []
    for path in sorted((ARTICLE / "supplement").glob("*.csv")):
        if path.stat().st_size <= MAX_DERIVED_BYTES:
            included.append(path)
        else:
            large_derived.append(path)
    unique = sorted(set(included), key=lambda path: str(path.relative_to(ROOT)))
    missing = [path for path in unique if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing reviewer files: {missing}")
    return unique, large_derived


def build_reviewer_archive(path: Path) -> None:
    included, large_derived = reviewer_files()
    manifest_rows = [
        (str(item.relative_to(ROOT)), item.stat().st_size, sha256(item)) for item in included
    ]
    readme_lines = [
        "# Peer-review code and derived summary tables",
        "",
        "This archive contains the analysis and figure-generation code, the exact public-input manifest, "
        "environment specification, manuscript sources, validation scripts and compact derived tables needed "
        "to audit the reported summary statistics.",
        "",
        "The third-party J-SHIS/NIED archives are not redistributed. Download locations, byte sizes and SHA-256 "
        "hashes are listed in public_inputs_manifest.tsv. Run work/verify_public_inputs.py before reproducing the analyses.",
        "",
        "Large station-level and spectrum-level derived matrices are excluded from this upload-size-controlled "
        "archive. They are deterministic outputs of the included scripts and will be included in the public "
        "repository/DOI release before publication.",
        "",
        "Excluded large derived files:",
    ]
    readme_lines.extend(
        f"- {item.relative_to(ROOT)} ({item.stat().st_size:,} bytes)" for item in large_derived
    )
    readme_lines.extend(
        [
            "",
            "Validation command:",
            "",
            "```bash",
            "python work/validate_compact_peer_review_archive.py",
            "```",
            "",
        ]
    )

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        archive.writestr("README_REVIEWERS.md", "\n".join(readme_lines))
        manifest = "path\tbytes\tsha256\n" + "\n".join(
            f"{relative}\t{size}\t{digest}" for relative, size, digest in manifest_rows
        ) + "\n"
        archive.writestr("MANIFEST_SHA256.tsv", manifest)
        for item in included:
            archive.write(item, arcname=str(item.relative_to(ROOT)))


def write_upload_manifest(files: list[tuple[str, Path]]) -> None:
    manifest_path = OUTPUT / "UPLOAD_MANIFEST.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["role", "file", "bytes", "sha256"])
        for role, path in files:
            writer.writerow([role, path.name, path.stat().st_size, sha256(path)])


def run() -> None:
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    OUTPUT.mkdir(parents=True)

    manuscript = OUTPUT / f"{PREFIX}_manuscript.pdf"
    supplement = OUTPUT / f"{PREFIX}_supplementary_information.pdf"
    cover = OUTPUT / f"{PREFIX}_cover_letter.txt"
    reviewer = OUTPUT / f"{PREFIX}_reviewer_code_and_derived_tables.zip"
    shutil.copy2(ARTICLE / "main.pdf", manuscript)
    shutil.copy2(ARTICLE / "supplementary_information.pdf", supplement)
    cover.write_text((ARTICLE / "cover_letter_cee.md").read_text(encoding="utf-8"), encoding="utf-8")
    build_reviewer_archive(reviewer)

    files = [
        ("Manuscript", manuscript),
        ("Supplementary Information", supplement),
        ("Cover letter", cover),
        ("Peer-review code and derived tables", reviewer),
    ]
    write_upload_manifest(files)
    instructions = [
        "CEE initial-submission files",
        "",
        "Upload the manuscript PDF as the main manuscript and the Supplementary Information PDF separately.",
        "Paste or upload the cover letter in the journal system.",
        "Upload the ZIP as a supplementary/related peer-review file if the portal permits; otherwise provide it when requested by the editor.",
        "UPLOAD_MANIFEST.tsv is a local integrity record and does not need to be uploaded.",
        "",
        "The GitHub repository is currently private. Make the reviewed release public or deposit it with a DOI before publication, then update the manuscript links if the final DOI changes.",
        "",
    ]
    (OUTPUT / "README_UPLOAD.txt").write_text("\n".join(instructions), encoding="utf-8")
    for role, path in files:
        print(f"{role}: {path.name} ({path.stat().st_size:,} bytes)")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    run()
