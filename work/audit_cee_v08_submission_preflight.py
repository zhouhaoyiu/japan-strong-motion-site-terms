#!/usr/bin/env python3
"""Preflight audit for the CEE v0.8 manuscript package."""

from __future__ import annotations

import csv
import re
import subprocess
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
MAIN = ARTICLE / "main.tex"
SUPP = ARTICLE / "supplementary_information.tex"
ZIP = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article.zip"
SHA = ROOT / "outputs" / "cee_submission_package_sha256.txt"
OUT_MD = ROOT / "outputs" / "cee_v08_submission_preflight_audit.md"
OUT_CSV = ROOT / "outputs" / "cee_v08_submission_preflight_audit.csv"


def tex_to_words(text: str) -> list[str]:
    text = text.replace(r"\%", " percent ")
    text = re.sub(r"%.*", " ", text)
    text = re.sub(r"\\url\{[^{}]*\}", " ", text)
    text = re.sub(r"\\path\{[^{}]*\}", " ", text)
    text = re.sub(r"\\textsuperscript\{[^{}]*\}", " ", text)
    text = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?", " ", text)
    text = re.sub(r"[{}$^_&~]", " ", text)
    return re.findall(r"[A-Za-z0-9]+(?:[-.][A-Za-z0-9]+)*", text)


def between(text: str, start: str, end: str) -> str:
    s = text.index(start)
    e = text.index(end, s)
    return text[s:e]


def section(text: str, name: str) -> str:
    pattern = re.compile(rf"\\section\*?\{{{re.escape(name)}\}}")
    match = pattern.search(text)
    if not match:
        return ""
    next_match = re.search(r"\\section\*?\{", text[match.end():])
    if not next_match:
        return text[match.end():]
    return text[match.end(): match.end() + next_match.start()]


def pdf_pages(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        result = subprocess.run(["pdfinfo", str(path)], check=True, capture_output=True, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    return None


def checksum_ok() -> bool:
    if not SHA.exists() or not ZIP.exists():
        return False
    try:
        subprocess.run(
            ["shasum", "-a", "256", "-c", str(SHA.relative_to(ROOT))],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError:
        return False
    return True


def status(ok: bool, severity: str = "must") -> str:
    if ok:
        return "pass"
    return "open_" + severity


def main() -> None:
    text = MAIN.read_text(encoding="utf-8")
    supp = SUPP.read_text(encoding="utf-8")

    title_match = re.search(r"\\LARGE\\bfseries\s+(.+?)\\par", text, re.S)
    title = re.sub(r"\s+", " ", title_match.group(1)).strip() if title_match else ""
    title_words = tex_to_words(title)

    abstract = between(text, r"\noindent\textbf{Abstract}", r"\noindent\textbf{Keywords:}")
    abstract = abstract.replace(r"\noindent\textbf{Abstract}", " ")
    abstract_words = tex_to_words(abstract)

    before_methods = text.split(r"\section{Data and methods}", 1)[0]
    main_body = text.split(r"\section*{Data availability}", 1)[0]
    body_words = tex_to_words(main_body)

    refs = re.findall(r"\\item\s+", section(text, "References"))
    main_figs = re.findall(r"\\begin\{figure\}", text)
    supp_figs = re.findall(r"\\begin\{figure\}", supp)

    data_avail = section(text, "Data availability")
    code_avail = section(text, "Code availability")
    ai_use = section(text, "Artificial intelligence use")
    author_contributions = section(text, "Author contributions")
    funding = section(text, "Funding")
    competing_interests = section(text, "Competing interests")
    acknowledgements = section(text, "Acknowledgements")
    author_block_ok = all(
        token in text
        for token in [
            "Haoyu Zhou",
            "Qiang Ma",
            "Institute of Engineering Mechanics",
            "maqiang@iem.ac.cn",
            "0009-0003-8817-1209",
            "0000-0002-9768-5223",
        ]
    )
    github_ok = "https://github.com/zhouhaoyiu/japan-strong-motion-site-terms" in text

    placeholder_patterns = {
        "final_archive_link_placeholder": r"final public archive citation|approved GitHub or Zenodo archive citation|will be added|before journal submission",
        "generic_placeholder_token": r"TBD|TODO|PLACEHOLDER|To be confirmed|to be confirmed",
    }
    open_placeholders = []
    for label, pattern in placeholder_patterns.items():
        if re.search(pattern, text, re.I):
            open_placeholders.append(label)

    phrase_patterns = {
        "rather_than": r"\brather than\b",
        "however": r"\bhowever\b",
        "therefore": r"\btherefore\b",
        "not_only": r"\bnot only\b",
        "not_x_but_y_candidate": r"\bnot\b.{0,80}\bbut\b",
    }
    phrase_hits = []
    for label, pattern in phrase_patterns.items():
        if re.search(pattern, text, re.I | re.S):
            phrase_hits.append(label)

    rows = [
        {
            "gate": "title_length",
            "status": status(0 < len(title_words) <= 15, "revise"),
            "evidence": f"{len(title_words)} words",
            "action": "Keep title concise for CEE editorial fit.",
        },
        {
            "gate": "abstract_length",
            "status": status(len(abstract_words) <= 150, "must"),
            "evidence": f"{len(abstract_words)} words",
            "action": "CEE Article abstracts should be 150 words or fewer.",
        },
        {
            "gate": "article_word_count_tracking",
            "status": status(len(body_words) <= 6500, "revise"),
            "evidence": f"{len(body_words)} words before availability/references",
            "action": "If target system enforces a strict count, reduce results/discussion density.",
        },
        {
            "gate": "data_availability",
            "status": status(bool(data_avail.strip()), "must"),
            "evidence": "present" if data_avail.strip() else "missing",
            "action": "Keep raw-data redistribution boundary explicit.",
        },
        {
            "gate": "code_availability",
            "status": status(bool(code_avail.strip()) and github_ok, "must"),
            "evidence": "present with GitHub URL" if code_avail.strip() and github_ok else "missing GitHub URL",
            "action": "Keep code repository URL visible in Code availability.",
        },
        {
            "gate": "author_metadata",
            "status": status(author_block_ok, "must"),
            "evidence": "present" if author_block_ok else "missing or incomplete",
            "action": "Confirm author names, affiliation, correspondence, and ORCID identifiers.",
        },
        {
            "gate": "ai_use_statement",
            "status": status(bool(ai_use.strip()), "must"),
            "evidence": "present" if ai_use.strip() else "missing",
            "action": "Keep journal-compliant AI-use disclosure.",
        },
        {
            "gate": "author_contributions",
            "status": status(bool(author_contributions.strip()), "must"),
            "evidence": "present" if author_contributions.strip() else "missing",
            "action": "Keep author contribution statement in the manuscript.",
        },
        {
            "gate": "funding_statement",
            "status": status(bool(funding.strip()) and "no specific grant" in funding.lower(), "must"),
            "evidence": "no specific funding declared" if funding.strip() else "missing",
            "action": "Keep no-funding statement unless coauthors add a grant.",
        },
        {
            "gate": "competing_interests",
            "status": status(bool(competing_interests.strip()) and "no competing interests" in competing_interests.lower(), "must"),
            "evidence": "none declared" if competing_interests.strip() else "missing",
            "action": "Keep competing-interest declaration.",
        },
        {
            "gate": "acknowledgements",
            "status": status(bool(acknowledgements.strip()), "must"),
            "evidence": "present" if acknowledgements.strip() else "missing",
            "action": "Keep data-provider acknowledgement details.",
        },
        {
            "gate": "references_count",
            "status": status(len(refs) >= 10, "revise"),
            "evidence": f"{len(refs)} references",
            "action": "Final reference-manager export remains required.",
        },
        {
            "gate": "main_and_supp_figures",
            "status": status(len(main_figs) >= 4 and len(supp_figs) >= 8, "revise"),
            "evidence": f"{len(main_figs)} main figures; {len(supp_figs)} supplementary figures",
            "action": "Final figure portal check remains required.",
        },
        {
            "gate": "placeholder_scan",
            "status": status(not open_placeholders, "must"),
            "evidence": "; ".join(open_placeholders) if open_placeholders else "none",
            "action": "Resolve final archive citation and any metadata placeholders before submission.",
        },
        {
            "gate": "style_guardrail_scan",
            "status": status(not phrase_hits, "revise"),
            "evidence": "; ".join(phrase_hits) if phrase_hits else "none",
            "action": "Revise only if a hit is real prose rather than necessary technical contrast.",
        },
        {
            "gate": "pdf_outputs",
            "status": status(pdf_pages(ARTICLE / "main.pdf") is not None and pdf_pages(ARTICLE / "supplementary_information.pdf") is not None, "must"),
            "evidence": f"main={pdf_pages(ARTICLE / 'main.pdf')} pages; supplement={pdf_pages(ARTICLE / 'supplementary_information.pdf')} pages",
            "action": "Recompile after final text edits.",
        },
        {
            "gate": "zip_checksum",
            "status": status(checksum_ok(), "must"),
            "evidence": "ok" if checksum_ok() else "outdated or missing",
            "action": "Rebuild zip and checksum after final package edits.",
        },
    ]

    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gate", "status", "evidence", "action"])
        writer.writeheader()
        writer.writerows(rows)

    open_rows = [row for row in rows if row["status"] != "pass"]
    lines = [
        "# CEE v0.8 submission preflight audit",
        "",
        f"Generated: {date.today().isoformat()}.",
        "",
        "This audit tracks CEE submission-readiness gates for the current LaTeX package. It is an internal QA file.",
        "",
        "## Summary",
        "",
        f"- Gates checked: {len(rows)}.",
        f"- Passing gates: {len(rows) - len(open_rows)}.",
        f"- Open gates: {len(open_rows)}.",
        "",
        "## Open gates",
        "",
    ]
    if open_rows:
        for row in open_rows:
            lines.append(f"- **{row['gate']}** ({row['status']}): {row['evidence']}. Action: {row['action']}")
    else:
        lines.append("- None.")
    lines.extend(["", "## Full gate table", ""])
    for row in rows:
        lines.append(f"- `{row['gate']}`: {row['status']} ({row['evidence']})")
    lines.append("")
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    print(f"wrote {OUT_CSV.relative_to(ROOT)}")
    for row in open_rows:
        print(f"open: {row['gate']} -> {row['evidence']}")


if __name__ == "__main__":
    main()
