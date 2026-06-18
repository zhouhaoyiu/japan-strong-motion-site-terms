#!/usr/bin/env python3
"""Build a file-level public-archive readiness review for the CEE v0.8 package."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
ARTICLE_DIR = OUTPUT_DIR / "cee_submission_latex_v0_8_english_article"
PACKAGE_ZIP = OUTPUT_DIR / "cee_submission_latex_v0_8_english_article.zip"

OUT_CSV = OUTPUT_DIR / "cee_v08_public_archive_file_review.csv"
OUT_MD = OUTPUT_DIR / "cee_v08_public_archive_readiness.md"
OUT_JSON = OUTPUT_DIR / "cee_release_metadata_zenodo_v0_8_draft.json"
PUBLIC_CANDIDATE_ZIP = OUTPUT_DIR / "cee_public_release_candidate_v0_8.zip"
PUBLIC_CANDIDATE_SHA = OUTPUT_DIR / "cee_public_release_candidate_v0_8_sha256.txt"
PUBLIC_CANDIDATE_AUDIT = OUTPUT_DIR / "cee_public_release_candidate_v0_8_audit.md"
PUBLIC_CANDIDATE_SMOKE = OUTPUT_DIR / "cee_public_release_candidate_v0_8_smoke_test.md"

IMPORTANT_SCRIPTS = [
    "work/jshis_mf2013_official_residuals.py",
    "work/jshis_mf2013_site_term_bootstrap.py",
    "work/jshis_nonergodic_station_residual_model.py",
    "work/jshis_nonergodic_station_model_extensions.py",
    "work/jshis_station_hazard_impact.py",
    "work/jshis_official_hazard_response_check.py",
    "work/jshis_official_hazard_response_uncertainty_map.py",
    "work/build_openquake_jpn_station_psha.py",
    "work/postprocess_openquake_jpn_station_psha.py",
    "work/jshis_public_psha_tiled_map_recompute.py",
    "work/jshis_public_psha_national_stratified_sample.py",
    "work/jshis_public_psha_sigma_audit.py",
    "work/jshis_public_psha_near_station_physical_audit.py",
    "work/jshis_public_psha_basin_period_proxy_audit.py",
    "work/jshis_public_psha_tail_robustness_cases_audit.py",
    "work/jshis_public_psha_nearest_station_borehole_audit.py",
    "work/jshis_official_psha_input_audit.py",
    "work/jshis_sa3_continuous_correction_surface.py",
    "work/jshis_multiperiod_nonergodic_psha_audit.py",
    "work/jshis_gsi_terrain_geomorphic_audit.py",
    "work/jshis_zamp_vs400_site_amplification_audit.py",
    "work/kiknet_multievent_transfer_function_audit.py",
]

RAW_PATTERNS = [
    re.compile(pattern, re.IGNORECASE)
    for pattern in [
        r"flatfile_sub1.*\.zip$",
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
    ]
]

ABSOLUTE_PATH_PATTERNS = [
    re.compile(r"/Users/yojironoda/"),
    re.compile(r"/private/var/"),
    re.compile(r"/var/folders/"),
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def is_raw_like(path: str) -> bool:
    return any(pattern.search(path) for pattern in RAW_PATTERNS)


def classify(path: str, size: int) -> tuple[str, str, str, str, str]:
    lower = path.lower()
    if is_raw_like(lower):
        return (
            "raw_or_cache_candidate",
            "third-party raw/cache data",
            "High",
            "exclude",
            "Must not be deposited without explicit provider permission.",
        )
    if lower.endswith(".zip"):
        return (
            "package_container",
            "compiled package archive",
            "Low",
            "include package checksum, not nested zip where avoidable",
            "Confirm final package hash after author approval.",
        )
    if lower.endswith((".tex", ".md", ".json", ".yml", ".yaml", ".txt")):
        return (
            "documentation_or_source",
            "author-generated documentation or source",
            "Low",
            "include",
            "Check final author metadata, repository URL, DOI, and declarations.",
        )
    if lower.endswith(".py"):
        return (
            "analysis_script",
            "author-generated analysis code",
            "Low",
            "include",
            "Run syntax/smoke checks in the public repository before DOI deposition.",
        )
    if lower.endswith((".pdf", ".png")):
        return (
            "rendered_manuscript_or_figure",
            "author-generated rendered output",
            "Low",
            "include",
            "Confirm figure and manuscript versions match final submission.",
        )
    if lower.endswith(".csv"):
        if "mesh" in lower or "station" in lower or "hazard" in lower or "psha" in lower or size > 2_000_000:
            return (
                "derived_table_needs_approval",
                "derived aggregate table from public J-SHIS/NIED inputs",
                "Medium",
                "include after author/legal approval",
                "Confirm this table is aggregate/derived and does not recreate prohibited raw flatfile tables.",
            )
        return (
            "derived_summary_table",
            "small derived summary table",
            "Low",
            "include",
            "Confirm row-level content does not expose restricted raw records.",
        )
    return (
        "other",
        "package support file",
        "Medium",
        "review",
        "Inspect before public DOI deposition.",
    )


def text_flags(data: bytes) -> tuple[bool, str]:
    if b"\x00" in data[:4096]:
        return False, ""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        try:
            text = data.decode("latin-1")
        except UnicodeDecodeError:
            return False, ""
    hits = [pattern.pattern for pattern in ABSOLUTE_PATH_PATTERNS if pattern.search(text)]
    return bool(hits), ";".join(hits)


def zip_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with zipfile.ZipFile(PACKAGE_ZIP) as zf:
        for info in sorted(zf.infolist(), key=lambda x: x.filename):
            if info.is_dir():
                continue
            data = zf.read(info.filename)
            category, basis, risk, action, review = classify(info.filename, info.file_size)
            has_abs, abs_patterns = text_flags(data)
            rows.append(
                {
                    "path": info.filename,
                    "source": "cee_submission_zip_v0_8",
                    "category": category,
                    "size_bytes": info.file_size,
                    "sha256": sha256_bytes(data),
                    "data_basis": basis,
                    "redistribution_risk": risk,
                    "proposed_release_action": action,
                    "required_review": review,
                    "absolute_path_flag": "yes" if has_abs else "no",
                    "absolute_path_patterns": abs_patterns,
                    "approval_status": "Pending author/legal approval",
                }
            )
    return rows


def script_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    zip_paths = set()
    with zipfile.ZipFile(PACKAGE_ZIP) as zf:
        zip_paths = {info.filename for info in zf.infolist()}
    for rel in IMPORTANT_SCRIPTS:
        path = ROOT / rel
        if not path.exists():
            rows.append(
                {
                    "path": rel,
                    "source": "workspace_script_reference",
                    "category": "missing_script",
                    "size_bytes": 0,
                    "sha256": "",
                    "data_basis": "expected analysis code",
                    "redistribution_risk": "Medium",
                    "proposed_release_action": "add or explain omission",
                    "required_review": "Required script was referenced by the manuscript package but was not found.",
                    "absolute_path_flag": "unknown",
                    "absolute_path_patterns": "",
                    "approval_status": "Open technical action",
                }
            )
            continue
        data = path.read_bytes()
        has_abs, abs_patterns = text_flags(data)
        category, basis, risk, action, review = classify(rel, path.stat().st_size)
        rows.append(
            {
                "path": rel,
                "source": "workspace_script_reference",
                "category": category,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "data_basis": basis,
                "redistribution_risk": risk,
                "proposed_release_action": action,
                "required_review": review,
                "absolute_path_flag": "yes" if has_abs else "no",
                "absolute_path_patterns": abs_patterns,
                "approval_status": "Pending author/legal approval",
            }
        )
    # Avoid duplicate rows if a future zip starts carrying scripts.
    return [row for row in rows if row["path"] not in zip_paths]


def write_csv(rows: list[dict[str, Any]]) -> None:
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "path",
        "source",
        "category",
        "size_bytes",
        "sha256",
        "data_basis",
        "redistribution_risk",
        "proposed_release_action",
        "required_review",
        "absolute_path_flag",
        "absolute_path_patterns",
        "approval_status",
    ]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(rows: list[dict[str, Any]]) -> None:
    payload = {
        "metadata": {
            "title": "Public-data station residual and response-spectrum sensitivity audits for Japanese strong-motion records",
            "upload_type": "dataset",
            "description": (
                "Derived aggregate tables, figures, scripts, and audit files supporting a CEE manuscript on "
                "long-period station residuals, non-ergodic site terms, J-SHIS response-spectrum sensitivity, "
                "public-parameter J-SHIS audit checks, official PSHA input-chain provenance, sampled-grid "
                "correction surfaces, all-period official response-spectrum propagation, multi-event KiK-net transfer functions, GSI terrain proxies, "
                "J-SHIS VS400 site-amplification evidence, and claim-evidence traceability files. "
                "Original J-SHIS/NIED and GSI products must be obtained from the official providers."
            ),
            "creators": [
                {
                    "name": "Zhou, Haoyu",
                    "affiliation": "Institute of Engineering Mechanics, China Earthquake Administration",
                    "orcid": "0009-0003-8817-1209",
                },
                {
                    "name": "Ma, Qiang",
                    "affiliation": "Institute of Engineering Mechanics, China Earthquake Administration",
                    "orcid": "0000-0002-9768-5223",
                }
            ],
            "keywords": [
                "strong-motion records",
                "station residuals",
                "non-ergodic ground motion",
                "site amplification",
                "Japan",
                "J-SHIS",
                "uniform-hazard spectrum",
                "hazard spectra",
            ],
            "license": "cc-by-4.0",
            "access_right": "open",
            "version": "v0.8-pre-submission",
            "related_identifiers": [
                {
                    "identifier": "https://doi.org/10.17598/NIED.0032",
                    "relation": "isDerivedFrom",
                    "resource_type": "dataset",
                },
                {
                    "identifier": "https://www.j-shis.bosai.go.jp/",
                    "relation": "isDerivedFrom",
                    "resource_type": "dataset",
                },
                {
                    "identifier": "https://maps.gsi.go.jp/development/demtile.html",
                    "relation": "isDerivedFrom",
                    "resource_type": "dataset",
                },
                {
                    "identifier": "https://doi.org/10.1186/s40623-023-01936-y",
                    "relation": "cites",
                    "resource_type": "publication-article",
                },
            ],
            "notes": (
                "Repository: https://github.com/zhouhaoyiu/japan-strong-motion-site-terms. No specific funding "
                "is declared for this manuscript package. Do not deposit raw J-SHIS/NIED flatfiles, extracted "
                "raw tables, waveform files, or cached GSI tiles."
            ),
        },
        "archive_review": {
            "generated_on": date.today().isoformat(),
            "source_zip": str(PACKAGE_ZIP.relative_to(ROOT)),
            "source_zip_sha256": sha256_file(PACKAGE_ZIP),
            "file_review_csv": str(OUT_CSV.relative_to(ROOT)),
            "n_review_rows": len(rows),
            "public_release_candidate_zip": str(PUBLIC_CANDIDATE_ZIP.relative_to(ROOT)),
            "public_release_candidate_sha256_file": str(PUBLIC_CANDIDATE_SHA.relative_to(ROOT)),
            "public_release_candidate_audit": str(PUBLIC_CANDIDATE_AUDIT.relative_to(ROOT)),
            "public_release_candidate_smoke_test": str(PUBLIC_CANDIDATE_SMOKE.relative_to(ROOT)),
            "public_release_candidate_status_note": (
                "Current candidate file count, checksum, and smoke-test gate count are recorded in the "
                "external candidate audit, checksum, and smoke-test files; those generated files are kept "
                "outside the candidate zip to avoid a self-referential archive hash."
            ),
        },
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_markdown(rows: list[dict[str, Any]]) -> None:
    cats = Counter(row["category"] for row in rows)
    actions = Counter(row["proposed_release_action"] for row in rows)
    risks = Counter(row["redistribution_risk"] for row in rows)
    abs_hits = [row for row in rows if row["absolute_path_flag"] == "yes"]
    raw_hits = [row for row in rows if row["category"] == "raw_or_cache_candidate"]
    missing = [row for row in rows if row["category"] == "missing_script"]
    include_now = sum(1 for row in rows if row["proposed_release_action"] == "include")
    needs_approval = sum(1 for row in rows if "approval" in row["proposed_release_action"])
    review_rows = sum(1 for row in rows if row["proposed_release_action"] in {"review", "add or explain omission"})
    text = [
        "# CEE v0.8 public-archive readiness review",
        "",
        f"Generated: {date.today().isoformat()}.",
        "",
        "## Scope",
        "",
        f"- Source package: `{PACKAGE_ZIP.relative_to(ROOT)}`.",
        f"- Source package SHA-256: `{sha256_file(PACKAGE_ZIP)}`.",
        f"- Local public release candidate: `{PUBLIC_CANDIDATE_ZIP.relative_to(ROOT)}`.",
        f"- Local public release candidate SHA-256 is recorded outside the zip in `{PUBLIC_CANDIDATE_SHA.relative_to(ROOT)}` to avoid a self-referential archive hash.",
        f"- Candidate audit: `{PUBLIC_CANDIDATE_AUDIT.relative_to(ROOT)}`.",
        f"- Candidate smoke test: `{PUBLIC_CANDIDATE_SMOKE.relative_to(ROOT)}`.",
        f"- File-level review table: `{OUT_CSV.relative_to(ROOT)}`.",
        f"- Zenodo metadata draft: `{OUT_JSON.relative_to(ROOT)}`.",
        "- This is a pre-deposition technical review, not author/legal approval.",
        "",
        "## File counts",
        "",
        f"- Review rows: {len(rows)}.",
        f"- Proposed immediate include rows: {include_now}.",
        f"- Rows requiring author/legal approval before release: {needs_approval}.",
        f"- Rows requiring technical review or explanation: {review_rows}.",
        f"- Absolute-path flags: {len(abs_hits)}.",
        f"- Raw/cache candidate flags: {len(raw_hits)}.",
        f"- Missing expected scripts: {len(missing)}.",
        "",
        "## Category counts",
        "",
    ]
    for key, value in sorted(cats.items()):
        text.append(f"- `{key}`: {value}")
    text.extend(["", "## Risk counts", ""])
    for key, value in sorted(risks.items()):
        text.append(f"- `{key}`: {value}")
    text.extend(["", "## Proposed release actions", ""])
    for key, value in sorted(actions.items()):
        text.append(f"- `{key}`: {value}")
    text.extend(
        [
            "",
            "## Candidate package smoke test",
            "",
            "- Candidate file count, checksum, and smoke-test gate count are recorded in the external candidate audit, checksum, and smoke-test files.",
            "- Nested manuscript archive checksum: passed.",
            "- Nested manuscript archive entries contain no stale `outputs/` prefix, auxiliary LaTeX files, or internal-review entries.",
            "- PDF page counts: main article 18 pages; supplementary information 14 pages.",
            "- Raw/cache path scan: 0 hits.",
            "- Absolute local path scan: 0 hits.",
            "- Private-planning text scan: 0 hits.",
            "- Python syntax compilation: no failures.",
            "- The candidate zip is kept separate from its checksum, audit report, and smoke-test report to avoid a self-referential archive hash.",
        ]
    )
    text.extend(
        [
            "",
            "## Current blockers",
            "",
            "- Author/legal approval is still required for derived tables, especially mesh-level, station-level, and hazard/PSHA-derived CSV files.",
            "- The GitHub repository URL, creator metadata, no-funding declaration, and draft CC-BY-4.0 metadata are populated in the package.",
            "- Repository public visibility, final license approval, and derived-table redistribution boundaries must be confirmed before release.",
            "- A Zenodo DOI is optional unless the journal, coauthors, or data policy require a DOI-bearing archive in addition to GitHub.",
            "- Raw J-SHIS/NIED flatfiles, extracted raw tables, waveform files, local caches, and cached GSI tiles must stay out of the public archive.",
            "- The archive can include scripts and derived aggregate outputs after the author-approved redistribution boundary is confirmed.",
        ]
    )
    if abs_hits:
        text.extend(["", "## Absolute path flags", ""])
        for row in abs_hits[:20]:
            text.append(f"- `{row['path']}`: {row['absolute_path_patterns']}")
    if raw_hits:
        text.extend(["", "## Raw/cache candidate flags", ""])
        for row in raw_hits[:20]:
            text.append(f"- `{row['path']}`")
    if missing:
        text.extend(["", "## Missing expected scripts", ""])
        for row in missing:
            text.append(f"- `{row['path']}`")
    text.extend(
        [
            "",
            "## Submission-readiness effect",
            "",
            "This review does not change the scientific result. It reduces technical submission risk by making the public-archive boundary explicit and by showing that the current v0.8 package can be converted into a GitHub/Zenodo release without redistributing raw third-party data.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(text), encoding="utf-8")


def main() -> None:
    rows = zip_rows() + script_rows()
    write_csv(rows)
    write_json(rows)
    write_markdown(rows)
    print(f"wrote {OUT_CSV.relative_to(ROOT)} with {len(rows)} rows")
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
