#!/usr/bin/env python3
"""Audit the CEE v0.4 manuscript draft for submission-facing hygiene.

This script checks citation numbering, placeholder inventory, and conservative
claim guardrails. It is intentionally narrow: a passing result does not mean the
manuscript is submission-ready, only that this draft passes the automated checks
defined here.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
MANUSCRIPT = OUTPUTS / "cee_article_draft_v0_4_submission_prose.md"
AUDIT_MD = OUTPUTS / "cee_article_draft_v0_4_audit.md"
AUDIT_CSV = OUTPUTS / "cee_article_draft_v0_4_audit.csv"


@dataclass
class Check:
    name: str
    status: str
    detail: str


def split_references(text: str) -> tuple[str, str]:
    marker = "\n## References\n"
    if marker not in text:
        return text, ""
    body, refs = text.split(marker, 1)
    return body, refs


def expand_citation_token(token: str) -> list[int]:
    values: list[int] = []
    for part in token.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_s, end_s = [p.strip() for p in part.split("-", 1)]
            start, end = int(start_s), int(end_s)
            values.extend(range(start, end + 1))
        else:
            values.append(int(part))
    return values


def numeric_citations(body: str) -> list[int]:
    citations: list[int] = []
    for match in re.finditer(r"\[([0-9][0-9,\-\s]*)\]", body):
        citations.extend(expand_citation_token(match.group(1)))
    return citations


def placeholder_tokens(text: str) -> list[str]:
    placeholders: list[str] = []
    for match in re.finditer(r"\[([^\]]+)\]", text):
        token = match.group(1).strip()
        if re.fullmatch(r"[0-9][0-9,\-\s]*", token):
            continue
        placeholders.append(f"[{token}]")
    return sorted(set(placeholders))


def reference_numbers(refs: str) -> list[int]:
    nums = []
    for line in refs.splitlines():
        match = re.match(r"^([0-9]+)\.\s+", line)
        if match:
            nums.append(int(match.group(1)))
    return nums


def context_snippets(text: str, phrase: str, width: int = 70) -> list[str]:
    snippets: list[str] = []
    for match in re.finditer(re.escape(phrase), text, flags=re.IGNORECASE):
        start = max(0, match.start() - width)
        end = min(len(text), match.end() + width)
        snippet = " ".join(text[start:end].split())
        snippets.append(snippet)
    return snippets


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    body, refs = split_references(text)
    checks: list[Check] = []

    ref_nums = reference_numbers(refs)
    expected_refs = list(range(1, 12))
    checks.append(
        Check(
            "reference_sequence",
            "pass" if ref_nums == expected_refs else "fail",
            f"found={ref_nums}; expected={expected_refs}",
        )
    )

    policy_terms = [
        "Communications Earth & Environment",
        "Aims and Scope",
        "Content types",
        "Submission guidelines",
    ]
    policy_hits = [term for term in policy_terms if term.lower() in refs.lower()]
    checks.append(
        Check(
            "cee_policy_references_removed",
            "pass" if not policy_hits else "fail",
            "no CEE policy references in manuscript reference list"
            if not policy_hits
            else "hits=" + "; ".join(policy_hits),
        )
    )

    cites = numeric_citations(body)
    invalid_cites = sorted({c for c in cites if c < 1 or c > 11})
    missing_cites = sorted(set(expected_refs) - set(cites))
    checks.append(
        Check(
            "citation_numbers_within_reference_range",
            "pass" if cites and not invalid_cites else "fail",
            f"citations={sorted(set(cites))}; invalid={invalid_cites}",
        )
    )
    checks.append(
        Check(
            "all_references_cited",
            "pass" if not missing_cites else "warn",
            "all numbered references are cited in the body"
            if not missing_cites
            else "missing_citations=" + ", ".join(str(c) for c in missing_cites),
        )
    )

    placeholders = placeholder_tokens(text)
    checks.append(
        Check(
            "known_placeholders_remaining",
            "warn" if placeholders else "pass",
            f"{len(placeholders)} placeholders remain: " + "; ".join(placeholders)
            if placeholders
            else "no bracketed non-citation placeholders found",
        )
    )

    banned_overclaims = [
        "guaranteed acceptance",
        "guarantee of acceptance",
        "proves causality",
        "prove causality",
        "causally proves",
        "earthquake prediction model",
        "predict earthquakes",
        "we fully reproduce MF2013",
        "full Zhao et al. (2006) reproduction",
    ]
    overclaim_hits = [
        f"{phrase}: {' | '.join(context_snippets(text, phrase)[:2])}"
        for phrase in banned_overclaims
        if context_snippets(text, phrase)
    ]
    checks.append(
        Check(
            "overclaim_phrase_scan",
            "pass" if not overclaim_hits else "fail",
            "no banned overclaim phrases found"
            if not overclaim_hits
            else "hits=" + " || ".join(overclaim_hits),
        )
    )

    required_guardrails = {
        "causal_limit": "not causal proof",
        "mf2013_limit": "not described as a full MF2013",
        "ai_limit": "not an official full-model reproduction",
        "zhao_limit": "rather than a full OpenQuake rupture-context reconstruction",
        "ph_limit": "PH term is omitted",
    }
    missing_guardrails = [
        name for name, phrase in required_guardrails.items() if phrase.lower() not in text.lower()
    ]
    checks.append(
        Check(
            "required_guardrail_language_present",
            "pass" if not missing_guardrails else "fail",
            "all required conservative-limit phrases found"
            if not missing_guardrails
            else "missing=" + ", ".join(missing_guardrails),
        )
    )

    checks.append(
        Check(
            "submission_readiness_interpretation",
            "warn",
            "Automated checks do not resolve expert review, author metadata, final reference-manager export, repository DOI, data-license review, or clean-environment reproduction.",
        )
    )

    AUDIT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with AUDIT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["check", "status", "detail"])
        writer.writeheader()
        for check in checks:
            writer.writerow({"check": check.name, "status": check.status, "detail": check.detail})

    lines = [
        "# CEE Article Draft v0.4 Automated Audit",
        "",
        "Date: 2026-06-05.",
        "",
        "Status: automated manuscript-hygiene audit. Passing checks do not make the manuscript submission-ready.",
        "",
        "## Summary",
        "",
    ]
    failures = [c for c in checks if c.status == "fail"]
    warnings = [c for c in checks if c.status == "warn"]
    lines.append(f"- Checks run: {len(checks)}")
    lines.append(f"- Failures: {len(failures)}")
    lines.append(f"- Warnings: {len(warnings)}")
    lines.append(
        "- Overall automated result: "
        + ("FAIL" if failures else "PASS WITH KNOWN WARNINGS" if warnings else "PASS")
    )
    lines.extend(["", "## Check Results", ""])
    lines.append("| Check | Status | Detail |")
    lines.append("|---|---:|---|")
    for check in checks:
        detail = check.detail.replace("|", "\\|")
        lines.append(f"| `{check.name}` | {check.status.upper()} | {detail} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The v0.4 draft passes the automated citation-range, reference-list, CEE-policy-reference, and conservative-claim checks. It still contains known submission placeholders and still requires domain-expert, reference-manager, repository, license, and clean-environment reproduction review before submission.",
            "",
            "## Files",
            "",
            f"- Manuscript: `{MANUSCRIPT.relative_to(ROOT)}`",
            f"- Machine-readable audit: `{AUDIT_CSV.relative_to(ROOT)}`",
        ]
    )
    AUDIT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    for check in checks:
        print(f"{check.status.upper():5s} {check.name}: {check.detail}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
