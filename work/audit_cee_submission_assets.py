#!/usr/bin/env python3
"""Audit CEE reference-import and submission-portal draft assets."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def first_existing(candidates: list[Path]) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


BIB = first_existing(
    [
        PROJECT_ROOT / "outputs" / "cee_reference_manager_import_v0_1.bib",
        PROJECT_ROOT / "docs" / "cee_reference_manager_import_v0_1.bib",
    ]
)
PORTAL = first_existing(
    [
        PROJECT_ROOT / "outputs" / "cee_submission_portal_fields_v0_1.md",
        PROJECT_ROOT / "docs" / "cee_submission_portal_fields_v0_1.md",
    ]
)
REPORT_DIR = PROJECT_ROOT / "outputs" if BIB.parent.name == "outputs" else PROJECT_ROOT / "qa"
REPORT = REPORT_DIR / "cee_submission_assets_audit.md"
CSV_REPORT = REPORT_DIR / "cee_submission_assets_audit.csv"


@dataclass
class Check:
    check: str
    status: str
    detail: str


def add(checks: list[Check], check: str, ok: bool, detail: str) -> None:
    checks.append(Check(check=check, status="PASS" if ok else "FAIL", detail=detail))


def display(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def bib_entries(text: str) -> list[tuple[str, str, str]]:
    pattern = re.compile(r"@(article|misc|inproceedings)\s*\{\s*([^,\s]+)\s*,", re.IGNORECASE)
    return [(match.group(1).lower(), match.group(2), match.group(0)) for match in pattern.finditer(text)]


def entry_body(text: str, key: str) -> str:
    start_match = re.search(r"@\w+\s*\{\s*" + re.escape(key) + r"\s*,", text)
    if not start_match:
        return ""
    start = start_match.start()
    next_match = re.search(r"\n@\w+\s*\{", text[start + 1 :])
    if next_match:
        return text[start : start + 1 + next_match.start()]
    return text[start:]


def audit_bib(checks: list[Check]) -> None:
    text = BIB.read_text(encoding="utf-8")
    entries = bib_entries(text)
    keys = [key for _, key, _ in entries]
    add(checks, "bib_entry_count", len(entries) == 14, f"entries={len(entries)}")
    add(checks, "bib_unique_keys", len(keys) == len(set(keys)), f"unique={len(set(keys))}; total={len(keys)}")

    required_doi_keys = [
        "NIED2025StrongMotionFlatFile2024",
        "MorikawaFujiwara2013MF2013",
        "StewartAfshariGoulet2017NonErgodicSiteResponse",
        "Lavrentiadis2023OverviewNonErgodicGMM",
        "Lacour2023EfficientNonErgodicLargeDatasets",
        "ParkerStewart2022ErgodicSiteResponseSubduction",
        "DeLaTorre2024WellingtonBasinResiduals",
        "Zhao2006AttenuationRelationsJapan",
        "Pagani2014OpenQuakeEngine",
    ]
    missing_doi = [key for key in required_doi_keys if "doi =" not in entry_body(text, key)]
    add(checks, "required_doi_fields", not missing_doi, "missing=" + ", ".join(missing_doi) if missing_doi else f"{len(required_doi_keys)} DOI fields present")

    ref2 = entry_body(text, "Morikawa2024FlatFileKNETKiKnet")
    ref2_ok = "verify full author list" in ref2.lower() and "no doi" in ref2.lower()
    add(checks, "ref2_manual_warning", ref2_ok, "Ref.2 warns about missing DOI and manual proceedings verification" if ref2_ok else "Ref.2 warning missing")

    policy_keys = ["CEE2026AimsScope", "CEE2026SubmissionGuidelines", "CEE2026ContentTypes"]
    policy_ok = all("not part of manuscript reference list" in entry_body(text, key).lower() for key in policy_keys)
    add(checks, "policy_reference_warning", policy_ok, "CEE policy entries marked as submission-planning sources")


def audit_portal(checks: list[Check]) -> None:
    text = PORTAL.read_text(encoding="utf-8")
    required_sections = [
        "## Manuscript Type",
        "## Cover Letter Text",
        "## Data Availability",
        "## Code Availability",
        "## AI/LLM Use Disclosure",
        "## Open Pre-Submission Fields",
        "## Integrity Note",
    ]
    missing = [section for section in required_sections if section not in text]
    add(checks, "portal_required_sections", not missing, "missing=" + ", ".join(missing) if missing else f"{len(required_sections)} sections present")

    guardrails = [
        "not as generic artificial-intelligence earthquake prediction",
        "not a full J-SHIS hazard-map reproduction",
        "does not claim causal proof from correlations",
        "PH Philippine Sea Plate term is explicitly omitted",
        "raw third-party data excluded",
    ]
    missing_guardrails = [phrase for phrase in guardrails if phrase not in text]
    add(checks, "portal_guardrail_language", not missing_guardrails, "missing=" + "; ".join(missing_guardrails) if missing_guardrails else f"{len(guardrails)} guardrails present")

    open_fields_ok = all(token in text for token in ["[Corresponding author name]", "[repository or DOI", "[license to be selected]", "Open |"])
    add(checks, "portal_open_placeholders_visible", open_fields_ok, "open placeholders remain visible for author/legal completion")


def write_reports(checks: list[Check]) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    with CSV_REPORT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"])
        writer.writeheader()
        for check in checks:
            writer.writerow(check.__dict__)

    failures = [check for check in checks if check.status != "PASS"]
    lines = [
        "# CEE Submission Assets Audit",
        "",
        "Date: 2026-06-05.",
        "",
        "Status: " + ("passed" if not failures else "failed"),
        "",
        "Scope: validates the draft BibTeX import scaffold and submission portal fields packet. Passing this audit does not replace final reference-manager verification, author sign-off, legal review, or journal portal checks.",
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    for check in checks:
        lines.append(f"| `{check.check}` | {check.status} | {check.detail} |")
    lines.extend(
        [
            "",
            "## Remaining Limits",
            "",
            "- Ref. 2 still requires final manual proceedings metadata and full author-list verification.",
            "- Author, funding, competing-interest, related-manuscript, repository/DOI, and license fields remain placeholders.",
            "- This audit checks consistency of draft assets only; it does not make the manuscript submission-ready.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    checks: list[Check] = []
    add(checks, "bib_file_exists", BIB.exists(), display(BIB))
    add(checks, "portal_file_exists", PORTAL.exists(), display(PORTAL))
    if BIB.exists():
        audit_bib(checks)
    if PORTAL.exists():
        audit_portal(checks)
    write_reports(checks)
    print(REPORT)
    if any(check.status != "PASS" for check in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
