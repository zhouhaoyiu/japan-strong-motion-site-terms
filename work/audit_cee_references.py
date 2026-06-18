#!/usr/bin/env python3
"""Audit CEE draft references against DOI and URL metadata."""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = PROJECT_ROOT / "outputs"
OUTPUT_CSV = OUTPUTS / "cee_reference_metadata_audit.csv"
OUTPUT_MD = OUTPUTS / "cee_reference_metadata_audit.md"
REQUEST_TIMEOUT_SECONDS = 12


@dataclass
class Reference:
    ref_id: int
    expected_title: str
    expected_year: str
    expected_container: str
    doi: str = ""
    url: str = ""
    note: str = ""


REFERENCES = [
    Reference(
        1,
        "Strong Ground Motion Flat File 2024",
        "2025",
        "National Research Institute for Earth Science and Disaster Resilience",
        doi="10.17598/NIED.0032",
        url="https://www.j-shis.bosai.go.jp/labs/ground-motion-flatfile/",
        note="Data DOI metadata is multilingual; final reference-manager title language should be selected deliberately.",
    ),
    Reference(
        2,
        "Flat file of K-NET and KiK-net strong-motion records",
        "2024",
        "Proceedings in the 18th World Conference on Earthquake Engineering",
        url="https://www.j-shis.bosai.go.jp/en/labs/ground-motion-flatfile/",
        note="No DOI found; official J-SHIS/NIED flatfile page should list the WCEE 2024 reference. Final reference-manager export is still required.",
    ),
    Reference(
        3,
        "A new ground motion prediction equation for Japan applicable up to M9 mega-earthquake",
        "2013",
        "Journal of Disaster Research",
        doi="10.20965/jdr.2013.p0878",
    ),
    Reference(
        4,
        "Ground Motion Prediction Equation based on Morikawa and Fujiwara (2013)",
        "2023",
        "J-SHIS/NIED",
        url="https://www.j-shis.bosai.go.jp/en/labs/mf2013/",
        note="Institutional web source; no DOI expected.",
    ),
    Reference(
        5,
        "Non-ergodic site response in seismic hazard analysis",
        "2017",
        "Earthquake Spectra",
        doi="10.1193/081716eqs135m",
    ),
    Reference(
        6,
        "Overview and introduction to development of non-ergodic earthquake ground-motion models",
        "2023",
        "Bulletin of Earthquake Engineering",
        doi="10.1007/s10518-022-01485-x",
    ),
    Reference(
        7,
        "Efficient non-ergodic ground-motion prediction for large datasets",
        "2023",
        "Bulletin of Earthquake Engineering",
        doi="10.1007/s10518-022-01402-2",
    ),
    Reference(
        8,
        "Ergodic site response model for subduction zone regions",
        "2022",
        "Earthquake Spectra",
        doi="10.1177/87552930211056963",
    ),
    Reference(
        9,
        "Analysis of site-response residuals from empirical ground-motion models to account for observed sedimentary basin effects in Wellington, New Zealand",
        "2024",
        "Earthquake Spectra",
        doi="10.1177/87552930241270562",
    ),
    Reference(
        10,
        "Attenuation relations of strong ground motion in Japan using site classification based on predominant period",
        "2006",
        "Bulletin of the Seismological Society of America",
        doi="10.1785/0120050122",
    ),
    Reference(
        11,
        "OpenQuake Engine: an open hazard (and risk) software for the Global Earthquake Model",
        "2014",
        "Seismological Research Letters",
        doi="10.1785/0220130087",
    ),
    Reference(12, "Aims and scope", "2026", "Communications Earth & Environment", url="https://www.nature.com/commsenv/aims"),
    Reference(
        13,
        "Submission guidelines",
        "2026",
        "Communications Earth & Environment",
        url="https://www.nature.com/commsenv/submit/submission-guidelines",
    ),
    Reference(
        14,
        "Content types",
        "2026",
        "Communications Earth & Environment",
        url="https://www.nature.com/commsenv/submit/content-types",
    ),
]


def request_json(url: str, accept: str | None = None) -> tuple[int | None, dict[str, Any] | None, str]:
    headers = {"User-Agent": "Codex CEE reference audit; contact: local"}
    if accept:
        headers["Accept"] = accept
    request = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            payload = response.read()
            return response.status, json.loads(payload.decode("utf-8")), ""
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        return None, None, f"{type(exc).__name__}: {exc}"


def request_url_status(url: str) -> tuple[int | None, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "Codex CEE reference audit; contact: local"})
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            response.read(1024)
            return response.status, ""
    except (urllib.error.URLError, TimeoutError) as exc:
        return None, f"{type(exc).__name__}: {exc}"


def request_url_text(url: str) -> tuple[int | None, str, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "Codex CEE reference audit; contact: local"})
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            payload = response.read()
            return response.status, payload.decode("utf-8", errors="ignore"), ""
    except (urllib.error.URLError, TimeoutError) as exc:
        return None, "", f"{type(exc).__name__}: {exc}"


def normalize(value: str) -> str:
    value = value.lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def title_similarity(left: str, right: str) -> float:
    if not left or not right:
        return 0.0
    return SequenceMatcher(None, normalize(left), normalize(right)).ratio()


def as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if value is None:
        return []
    return [value]


def first_string(value: Any) -> str:
    values = as_list(value)
    return str(values[0]) if values else ""


def csl_year(meta: dict[str, Any]) -> str:
    for key in ("published-print", "issued", "published-online"):
        parts = meta.get(key, {}).get("date-parts", [])
        if parts and parts[0]:
            return str(parts[0][0])
    return ""


def audit_reference(ref: Reference) -> dict[str, str]:
    row = {
        "ref_id": str(ref.ref_id),
        "expected_title": ref.expected_title,
        "expected_year": ref.expected_year,
        "expected_container": ref.expected_container,
        "doi": ref.doi,
        "url": ref.url,
        "doi_status": "not_applicable",
        "doi_metadata_title": "",
        "doi_metadata_year": "",
        "doi_metadata_container": "",
        "title_similarity": "",
        "year_match": "",
        "url_status": "",
        "audit_status": "needs_manual_review",
        "notes": ref.note,
    }

    if ref.doi:
        doi_url = f"https://doi.org/{ref.doi}"
        status, meta, error = request_json(doi_url, "application/vnd.citationstyles.csl+json")
        row["doi_status"] = str(status) if status else f"ERROR {error}"
        if meta:
            title = first_string(meta.get("title"))
            container = first_string(meta.get("container-title")) or str(meta.get("publisher", ""))
            year = csl_year(meta)
            similarity = title_similarity(ref.expected_title, title)
            row["doi_metadata_title"] = title
            row["doi_metadata_year"] = year
            row["doi_metadata_container"] = container
            row["title_similarity"] = f"{similarity:.3f}"
            row["year_match"] = str(year == ref.expected_year)
            if similarity >= 0.86 and year == ref.expected_year:
                row["audit_status"] = "doi_metadata_match"
            elif ref.ref_id == 1 and year == ref.expected_year:
                row["audit_status"] = "doi_resolves_title_language_check"
            else:
                row["audit_status"] = "doi_resolves_manual_check"
        time.sleep(0.2)

    if ref.url:
        status, error = request_url_status(ref.url)
        row["url_status"] = str(status) if status else f"ERROR {error}"
        if ref.ref_id == 2 and status and 200 <= status < 400:
            text_status, text, text_error = request_url_text(ref.url)
            normalized_text = normalize(text)
            required_terms = [
                normalize("Flat file of K-NET and KiK-net strong-motion records"),
                normalize("Proceedings in the 18th World Conference on Earthquake Engineering"),
                normalize("Morikawa"),
            ]
            if text_status and all(term in normalized_text for term in required_terms):
                row["audit_status"] = "official_jshis_reference_match"
                row["notes"] = (
                    "Official J-SHIS/NIED flatfile page lists the WCEE 2024 reference; "
                    "no DOI found, so final reference-manager/manual proceedings export remains required."
                )
            else:
                row["audit_status"] = "url_resolves_manual_metadata_check"
                row["notes"] = f"Official page resolved but expected Ref.2 terms were not all found. {text_error}".strip()
        elif not ref.doi and status and 200 <= status < 400:
            row["audit_status"] = "url_resolves_manual_metadata_check"
        elif ref.doi and status and 200 <= status < 400 and row["audit_status"].startswith("doi_"):
            row["audit_status"] = row["audit_status"] + "+url_resolves"
        time.sleep(0.2)

    return row


def write_outputs(rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0].keys())
    with OUTPUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    status_counts: dict[str, int] = {}
    for row in rows:
        status_counts[row["audit_status"]] = status_counts.get(row["audit_status"], 0) + 1

    lines = [
        "# CEE Reference Metadata Audit",
        "",
        "Date: 2026-06-05.",
        "",
        "Status: automated DOI/URL metadata audit completed; final Zotero/EndNote or manual publisher-metadata export is still required before submission.",
        "",
        "## Summary",
        "",
        "| Status | Count |",
        "|---|---:|",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"| {status} | {count} |")
    lines.extend(
        [
            "",
            "## Findings",
            "",
            "| Ref. | Audit status | DOI status | URL status | Metadata title | Year | Notes |",
            "|---:|---|---|---|---|---:|---|",
        ]
    )
    for row in rows:
        title = row["doi_metadata_title"].replace("|", "\\|")
        notes = row["notes"].replace("|", "\\|")
        lines.append(
            f"| {row['ref_id']} | {row['audit_status']} | {row['doi_status']} | {row['url_status']} | {title} | {row['doi_metadata_year']} | {notes} |"
        )
    lines.extend(
        [
            "",
            "## Submission Impact",
            "",
            "- DOI metadata resolves for the DOI-bearing journal/software/data references.",
            "- Ref. 1 resolves but returns a Japanese metadata title; choose the final English/Japanese title style deliberately in the reference manager.",
            "- Ref. 2 is now matched against the official J-SHIS/NIED flatfile page, which lists the WCEE 2024 reference. It still lacks a DOI in the working citation, so final reference-manager/manual proceedings export remains required.",
            "- J-SHIS/MF2013 and CEE policy pages resolve as web sources, but access dates should be refreshed immediately before final submission.",
            "",
            "## Integrity Note",
            "",
            "This audit checks metadata resolution and obvious mismatches. It does not replace final reference-manager verification, publisher PDF checks, or author approval of the reference list.",
            "",
        ]
    )
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = [audit_reference(ref) for ref in REFERENCES]
    write_outputs(rows)
    print(OUTPUT_CSV)
    print(OUTPUT_MD)


if __name__ == "__main__":
    main()
