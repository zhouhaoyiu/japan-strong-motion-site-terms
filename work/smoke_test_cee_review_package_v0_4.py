#!/usr/bin/env python3
"""Smoke-test the CEE v0.4 internal review package.

This checks package integrity and review evidence consistency without reading
raw J-SHIS/NIED flatfiles or local waveform/cache data.
"""

from __future__ import annotations

import hashlib
import py_compile
import shutil
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ZIP = PROJECT_ROOT / "outputs" / "cee_review_package_v0_4_internal.zip"
DEFAULT_EXTRACT_DIR = PROJECT_ROOT / "work" / "cee_review_package_v0_4_internal_smoke"
DEFAULT_REPORT = PROJECT_ROOT / "outputs" / "cee_review_package_v0_4_internal_smoke_test_report.md"

PACKAGE_DIRNAME = "cee_review_package_v0_4_internal"
FORBIDDEN_NAME_TOKENS = [
    "flatfile",
    "smrec_schema",
    "site_schema.tsv",
    "source_schema.tsv",
    ".mseed",
    ".BSON",
    "seisbench",
    "external_data",
    "__pycache__",
    ".pyc",
]
FORBIDDEN_TEXT_TOKENS = ["/" + "Users/", "/" + "private/", "/" + "var/folders"]


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


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def png_dimensions(path: Path) -> tuple[int, int]:
    header = path.read_bytes()[:24]
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path} is not a PNG")
    return struct.unpack(">II", header[16:24])


def add(checks: list[Check], name: str, ok: bool, detail: str) -> None:
    checks.append(Check(name=name, status="PASS" if ok else "FAIL", detail=detail))


def expected_file_checks(pkg: Path, checks: list[Check]) -> None:
    expected = [
        "README.md",
        "PACKAGE_MANIFEST.txt",
        "release/SHA256SUMS.txt",
        "release/environment.yml",
        "manuscript/cee_article_draft_v0_4_submission_prose.md",
        "manuscript/cee_combined_manuscript_v0_4_review.docx",
        "manuscript/cee_combined_manuscript_v0_4_review.pdf",
        "docs/cee_article_draft_v0_4_audit.md",
        "docs/cee_article_draft_v0_4_audit.csv",
        "docs/cee_ref2_manual_source_check.md",
        "docs/cee_domain_expert_review_packet.md",
        "docs/cee_reference_manager_verification_worksheet.md",
        "docs/cee_reference_manager_verification_checks.csv",
        "docs/cee_public_release_file_review.md",
        "docs/cee_public_release_file_review.csv",
        "docs/cee_scientific_review_worksheet.md",
        "docs/cee_scientific_review_checks.csv",
        "docs/cee_placeholder_resolution_worksheet.md",
        "docs/cee_placeholder_inventory.csv",
        "docs/cee_figure_production_audit.md",
        "docs/cee_figure_production_audit.csv",
        "docs/cee_staged_ai_reproduction_report.md",
        "docs/cee_staged_ai_reproduction_comparison.csv",
        "docs/cee_staged_zhao_reproduction_report.md",
        "docs/cee_staged_zhao_reproduction_comparison.csv",
        "docs/cee_fresh_clone_public_reproduction_report.md",
        "docs/cee_submission_action_packet.md",
        "docs/cee_submission_action_items.csv",
        "docs/cee_submission_blocker_dashboard.md",
        "docs/cee_submission_blocker_dashboard.csv",
        "docs/cee_submission_blocker_dashboard_audit.md",
        "docs/cee_submission_blocker_dashboard_audit.csv",
        "docs/cee_quantitative_claim_consistency_audit.md",
        "docs/cee_quantitative_claim_consistency_audit.csv",
        "docs/cee_submission_signoff_forms.md",
        "docs/cee_submission_signoff_register.csv",
        "docs/cee_submission_request_templates.md",
        "docs/cee_submission_request_tracker.csv",
        "docs/cee_submission_metadata_collection_form.md",
        "docs/cee_submission_metadata_required_fields.csv",
        "scripts/audit_cee_article_v0_4.py",
        "scripts/build_cee_review_docx_v0_4.py",
        "scripts/audit_cee_figure_production.py",
        "qa/docx_render_qa_report.md",
        "qa/docx_render/contact_sheet.png",
        "qa/package_integrity_qa_report.md",
        "tables/cee_claims_evidence_table.csv",
        "tables/jshis_mf2013_site_term_bootstrap_ci.csv",
        "tables/jshis_zhao2006_external_gmpe_station_site_correlations.csv",
        "tables/cee_fig5_zhao2006_external_gmpe_sensitivity_source.csv",
        "figures/cee_fig1_workflow_data_evidence.png",
        "figures/cee_fig2_mf2013_site_term_gain_ci.png",
        "figures/cee_fig3_residual_site_correlation_shift.png",
        "figures/cee_fig4_subgroup_robustness.png",
        "figures/cee_fig5_zhao2006_external_gmpe_sensitivity.png",
    ]
    missing = [item for item in expected if not (pkg / item).exists()]
    add(checks, "expected files", not missing, "missing: " + ", ".join(missing) if missing else f"{len(expected)} required files present")


def forbidden_scan(pkg: Path, checks: list[Check]) -> None:
    files = [path for path in pkg.rglob("*") if path.is_file()]
    name_hits = [
        str(path.relative_to(pkg))
        for path in files
        if any(token.lower() in str(path.relative_to(pkg)).lower() for token in FORBIDDEN_NAME_TOKENS)
    ]
    add(checks, "forbidden filenames", not name_hits, "hits: " + ", ".join(name_hits[:10]) if name_hits else "none")

    text_hits: list[str] = []
    for path in files:
        if path.suffix.lower() in {".md", ".py", ".txt", ".yml", ".yaml", ".csv"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if any(token in text for token in FORBIDDEN_TEXT_TOKENS):
                text_hits.append(str(path.relative_to(pkg)))
    add(checks, "local absolute paths", not text_hits, "hits: " + ", ".join(text_hits[:10]) if text_hits else "none")


def script_compile_checks(pkg: Path, checks: list[Check]) -> None:
    failures: list[str] = []
    scripts = sorted((pkg / "scripts").glob("*.py"))
    compile_dir = pkg.parent / "_py_compile"
    compile_dir.mkdir(parents=True, exist_ok=True)
    for script in scripts:
        try:
            py_compile.compile(str(script), cfile=str(compile_dir / f"{script.stem}.pyc"), doraise=True)
        except py_compile.PyCompileError as exc:
            failures.append(f"{script.name}: {exc.msg}")
    add(checks, "script compilation", not failures, "failures: " + "; ".join(failures) if failures else f"{len(scripts)} scripts compile")


def table_checks(pkg: Path, checks: list[Check]) -> None:
    tables = pkg / "tables"
    expectations = {
        "cee_claims_evidence_table.csv": (10, 7),
        "jshis_mf2013_site_term_bootstrap_ci.csv": (28, 19),
        "jshis_zhao2006_external_gmpe_station_site_correlations.csv": (32, 8),
        "cee_fig5_zhao2006_external_gmpe_sensitivity_source.csv": (31, 7),
        "jshis_zhao2006_external_gmpe_coverage_summary.csv": (4, 4),
    }
    failures: list[str] = []
    for name, expected_shape in expectations.items():
        df = pd.read_csv(tables / name)
        if df.shape != expected_shape:
            failures.append(f"{name}: expected {expected_shape}, got {df.shape}")
    add(checks, "table shapes", not failures, "failures: " + "; ".join(failures) if failures else f"{len(expectations)} expected table shapes match")

    boot = pd.read_csv(tables / "jshis_mf2013_site_term_bootstrap_ci.csv")
    all_event = boot[(boot["cluster_level"].eq("event")) & (boot["subgroup"].eq("all"))].copy()
    sa3 = all_event[all_event["target_label"].eq("SA(3.0s) RotD50")]["percent_mae_reduction"].iloc[0]
    sa1 = all_event[all_event["target_label"].eq("SA(1.0s) RotD50")]["percent_mae_reduction"].iloc[0]
    add(
        checks,
        "MF2013 key reductions",
        abs(sa1 - 16.067154) < 0.001 and abs(sa3 - 21.169213) < 0.001,
        f"SA1 event/all={sa1:.6f}; SA3 event/all={sa3:.6f}",
    )

    zhao = pd.read_csv(tables / "jshis_zhao2006_external_gmpe_station_site_correlations.csv")
    zhao_sa3_d1400 = zhao[(zhao["target_label"].eq("SA(3.0s) RotD50")) & (zhao["feature"].eq("log10(D1400)"))][
        "spearman_rho"
    ].iloc[0]
    add(
        checks,
        "Zhao D1400 key signal",
        abs(zhao_sa3_d1400 - 0.512581) < 0.001,
        f"SA3/D1400 rho={zhao_sa3_d1400:.6f}",
    )


def figure_and_document_checks(pkg: Path, checks: list[Check]) -> None:
    figures = pkg / "figures"
    expected_png_sizes = {
        "cee_fig1_workflow_data_evidence.png": (3150, 1680),
        "cee_fig2_mf2013_site_term_gain_ci.png": (2640, 1560),
        "cee_fig3_residual_site_correlation_shift.png": (3000, 1560),
        "cee_fig4_subgroup_robustness.png": (1920, 1560),
        "cee_fig5_zhao2006_external_gmpe_sensitivity.png": (3960, 1620),
    }
    failures: list[str] = []
    for name, expected_size in expected_png_sizes.items():
        size = png_dimensions(figures / name)
        if size != expected_size:
            failures.append(f"{name}: expected {expected_size}, got {size}")
        pdf = figures / name.replace(".png", ".pdf")
        if not pdf.read_bytes().startswith(b"%PDF"):
            failures.append(f"{pdf.name}: missing PDF header")
    add(checks, "figure files", not failures, "failures: " + "; ".join(failures) if failures else "5 PNG/PDF figure pairs verified")

    docx = pkg / "manuscript" / "cee_combined_manuscript_v0_4_review.docx"
    pdf = pkg / "manuscript" / "cee_combined_manuscript_v0_4_review.pdf"
    page_pngs = sorted((pkg / "qa" / "docx_render").glob("page-*.png"))
    ok = zipfile.is_zipfile(docx) and pdf.read_bytes().startswith(b"%PDF") and len(page_pngs) == 16
    add(checks, "review document files", ok, f"DOCX zip={zipfile.is_zipfile(docx)}; PDF header={pdf.read_bytes().startswith(b'%PDF')}; QA pages={len(page_pngs)}")

    qa_report = (pkg / "qa" / "docx_render_qa_report.md").read_text(encoding="utf-8")
    numbering_fixed = "references 8-11 correctly" in qa_report
    add(checks, "v0.4 reference-numbering QA", numbering_fixed, "QA report records explicit references 8-11 fix" if numbering_fixed else "reference-numbering fix note missing")


def figure_production_audit_checks(pkg: Path, checks: list[Check]) -> None:
    report = pkg / "docs" / "cee_figure_production_audit.md"
    csv_report = pkg / "docs" / "cee_figure_production_audit.csv"
    if not report.exists() or not csv_report.exists():
        add(checks, "figure production audit", False, "report or CSV missing")
        return
    report_text = report.read_text(encoding="utf-8")
    audit = pd.read_csv(csv_report)
    statuses = sorted(set(audit["status"]))
    ok = "Status: passed" in report_text and statuses == ["PASS"] and len(audit) == 41
    detail = f"report status={'passed' if 'Status: passed' in report_text else 'not passed'}; audit rows={len(audit)}; statuses={statuses}"
    add(checks, "figure production audit", ok, detail)


def staged_ai_reproduction_checks(pkg: Path, checks: list[Check]) -> None:
    report = pkg / "docs" / "cee_staged_ai_reproduction_report.md"
    comparison = pkg / "docs" / "cee_staged_ai_reproduction_comparison.csv"
    if not report.exists() or not comparison.exists():
        add(checks, "staged AI reproduction audit", False, "report or comparison CSV missing")
        return

    report_text = report.read_text(encoding="utf-8")
    audit = pd.read_csv(comparison)
    statuses = sorted(set(audit["status"]))
    ok = "Status: passed." in report_text and statuses == ["PASS"] and len(audit) == 2
    detail = f"report status={'passed' if 'Status: passed.' in report_text else 'not passed'}; audit rows={len(audit)}; statuses={statuses}"
    add(checks, "staged AI reproduction audit", ok, detail)


def staged_zhao_reproduction_checks(pkg: Path, checks: list[Check]) -> None:
    report = pkg / "docs" / "cee_staged_zhao_reproduction_report.md"
    comparison = pkg / "docs" / "cee_staged_zhao_reproduction_comparison.csv"
    if not report.exists() or not comparison.exists():
        add(checks, "staged Zhao reproduction audit", False, "report or comparison CSV missing")
        return

    report_text = report.read_text(encoding="utf-8")
    audit = pd.read_csv(comparison)
    statuses = sorted(set(audit["status"]))
    ok = "Status: passed." in report_text and statuses == ["PASS"] and len(audit) == 12
    detail = f"report status={'passed' if 'Status: passed.' in report_text else 'not passed'}; audit rows={len(audit)}; statuses={statuses}"
    add(checks, "staged Zhao reproduction audit", ok, detail)


def fresh_clone_reproduction_checks(pkg: Path, checks: list[Check]) -> None:
    report = pkg / "docs" / "cee_fresh_clone_public_reproduction_report.md"
    if not report.exists():
        add(checks, "fresh clone public package QA", False, "fresh-clone report missing")
        return

    text = report.read_text(encoding="utf-8")
    required = [
        "Status: passed.",
        "clean local clone",
        "not proof of a GitHub release asset or DOI deposition",
        "zip member CRC: PASS",
    ]
    missing = [item for item in required if item not in text]
    add(
        checks,
        "fresh clone public package QA",
        not missing,
        "clean-clone public package reproduction documented" if not missing else "; ".join(missing),
    )


def submission_process_docs_checks(pkg: Path, checks: list[Check]) -> None:
    docs = pkg / "docs"
    csv_expectations = {
        "cee_submission_action_items.csv": ((12, 10), "status", ["Open"]),
        "cee_submission_blocker_dashboard.csv": ((12, 12), "status", ["Open"]),
        "cee_submission_blocker_dashboard_audit.csv": ((24, 7), "status", ["PASS"]),
        "cee_quantitative_claim_consistency_audit.csv": ((66, 9), "status", ["PASS"]),
        "cee_submission_signoff_register.csv": ((6, 8), "status", ["Pending"]),
        "cee_submission_request_tracker.csv": ((5, 7), "status", ["Not sent"]),
        "cee_submission_metadata_required_fields.csv": ((20, 8), "status", ["Open"]),
        "cee_reference_manager_verification_checks.csv": ((14, 7), "status", ["Open"]),
        "cee_public_release_file_review.csv": ((101, 7), "approval_status", ["Pending author/legal approval"]),
        "cee_scientific_review_checks.csv": ((18, 9), "status", ["Pending expert review"]),
        "cee_placeholder_inventory.csv": ((65, 10), "status", ["Open"]),
    }
    failures: list[str] = []
    for filename, (expected_shape, status_column, expected_statuses) in csv_expectations.items():
        df = pd.read_csv(docs / filename)
        statuses = sorted(set(df[status_column]))
        if df.shape != expected_shape or statuses != expected_statuses:
            failures.append(f"{filename}: shape={df.shape}, statuses={statuses}")
    add(
        checks,
        "submission process tables",
        not failures,
        "action/signoff/request/metadata schemas verified" if not failures else "; ".join(failures),
    )

    text_expectations = {
        "cee_submission_action_packet.md": "These percentages are planning estimates, not publication-probability claims.",
        "cee_submission_blocker_dashboard.md": [
            "This dashboard is not approval",
            "is not a publication-probability estimate",
            "Do not submit while any blocking gate remains open.",
        ],
        "cee_submission_blocker_dashboard_audit.md": [
            "This audit is not approval",
            "All dashboard consistency checks currently pass.",
        ],
        "cee_quantitative_claim_consistency_audit.md": [
            "This audit is not peer review",
            "All checked quantitative manuscript claims currently match their source evidence.",
        ],
        "cee_submission_signoff_forms.md": "A blank or partially filled form is not approval.",
        "cee_submission_request_templates.md": "Do not state that the manuscript is submission-ready.",
        "cee_submission_metadata_collection_form.md": "Blank fields are not approval, and placeholder text must not be copied into the submission portal.",
        "cee_reference_manager_verification_worksheet.md": "This is not a verified reference export.",
        "cee_public_release_file_review.md": "This is not author/legal approval",
        "cee_scientific_review_worksheet.md": "This is not expert approval.",
        "cee_domain_expert_review_packet.md": "This file is not a reviewer report",
        "cee_placeholder_resolution_worksheet.md": [
            "This is not placeholder resolution",
            "Do not submit while required placeholders remain.",
        ],
    }
    missing: list[str] = []
    for filename, required_texts in text_expectations.items():
        if isinstance(required_texts, str):
            required_texts = [required_texts]
        text = (docs / filename).read_text(encoding="utf-8")
        for required_text in required_texts:
            if required_text not in text:
                missing.append(f"{filename}: {required_text}")
    add(
        checks,
        "submission process guardrails",
        not missing,
        "readiness/signoff/request guardrails present" if not missing else "; ".join(missing),
    )


def v04_manuscript_audit_checks(pkg: Path, checks: list[Check]) -> None:
    audit = pd.read_csv(pkg / "docs" / "cee_article_draft_v0_4_audit.csv")
    failures = audit[audit["status"].eq("fail")]
    warnings = audit[audit["status"].eq("warn")]
    required = {
        "reference_sequence",
        "cee_policy_references_removed",
        "citation_numbers_within_reference_range",
        "all_references_cited",
        "known_placeholders_remaining",
        "overclaim_phrase_scan",
        "required_guardrail_language_present",
        "submission_readiness_interpretation",
    }
    found = set(audit["check"])
    missing = sorted(required - found)
    ok = failures.empty and not missing and len(warnings) == 2
    detail = (
        f"checks={len(audit)}; failures={len(failures)}; warnings={len(warnings)}"
        if ok
        else f"missing={missing}; failures={len(failures)}; warnings={len(warnings)}"
    )
    add(checks, "v0.4 manuscript audit", ok, detail)


def reference_audit_checks(pkg: Path, checks: list[Check]) -> None:
    ref_audit = pd.read_csv(pkg / "docs" / "cee_reference_metadata_audit.csv")
    ref2 = ref_audit[ref_audit["ref_id"].eq(2)]
    if ref2.empty:
        add(checks, "Ref.2 official source check", False, "Ref.2 row missing from reference audit")
        return
    status = ref2["audit_status"].iloc[0]
    note = ref2["notes"].iloc[0]
    ok = status == "official_jshis_reference_match" and "J-SHIS/NIED" in note
    detail = f"status={status}; note={note}"
    add(checks, "Ref.2 official source check", ok, detail)


def manifest_checksum_checks(pkg: Path, checks: list[Check]) -> None:
    manifest = (pkg / "PACKAGE_MANIFEST.txt").read_text(encoding="utf-8").splitlines()
    actual = sorted(str(path.relative_to(pkg)) for path in pkg.rglob("*") if path.is_file())
    add(checks, "manifest matches files", manifest == actual, f"manifest={len(manifest)} files; actual={len(actual)} files")

    checksum_lines = (pkg / "release" / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    checksum_map = {}
    for line in checksum_lines:
        digest, rel = line.split(maxsplit=1)
        checksum_map[rel] = digest
    failures: list[str] = []
    for rel, expected_digest in checksum_map.items():
        path = pkg / rel
        if not path.exists():
            failures.append(f"missing {rel}")
        elif sha256(path) != expected_digest:
            failures.append(f"checksum mismatch {rel}")
    add(checks, "per-file SHA256", not failures, "failures: " + "; ".join(failures[:10]) if failures else f"{len(checksum_map)} checksums verified")


def write_report(path: Path, zip_path: Path, checks: list[Check]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    failures = [check for check in checks if check.status != "PASS"]
    lines = [
        "# CEE v0.4 Internal Review Package Smoke Test",
        "",
        "Date: 2026-06-06.",
        "",
        f"Package zip: `{display_path(zip_path)}`.",
        f"Package SHA256: `{sha256(zip_path)}`.",
        "",
        "Status: " + ("passed" if not failures else "failed"),
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    for check in checks:
        lines.append(f"| {check.name} | {check.status} | {check.detail} |")
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "This smoke test verifies internal review package integrity and derived-output consistency. It does not rerun the full raw-data MF2013 or Zhao workflows because the raw J-SHIS/NIED flatfile is intentionally excluded from the package.",
            "",
            "## Integrity Note",
            "",
            "Passing this smoke test does not make the manuscript submission-ready. Domain-expert review, reference-manager verification, public repository/DOI setup, and author/legal review remain required.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    zip_path = DEFAULT_ZIP
    extract_root = DEFAULT_EXTRACT_DIR
    report_path = DEFAULT_REPORT

    if extract_root.exists():
        shutil.rmtree(extract_root)
    extract_root.mkdir(parents=True)

    checks: list[Check] = []
    add(checks, "zip exists", zip_path.exists(), display_path(zip_path))
    add(checks, "zip integrity", zipfile.is_zipfile(zip_path), "zipfile header valid")

    with zipfile.ZipFile(zip_path) as archive:
        bad_member = archive.testzip()
        add(checks, "zip member CRC", bad_member is None, "all members OK" if bad_member is None else f"bad member: {bad_member}")
        archive.extractall(extract_root)

    pkg = extract_root / PACKAGE_DIRNAME
    add(checks, "extracted package directory", pkg.exists(), display_path(pkg))

    expected_file_checks(pkg, checks)
    forbidden_scan(pkg, checks)
    script_compile_checks(pkg, checks)
    table_checks(pkg, checks)
    figure_and_document_checks(pkg, checks)
    figure_production_audit_checks(pkg, checks)
    staged_ai_reproduction_checks(pkg, checks)
    staged_zhao_reproduction_checks(pkg, checks)
    fresh_clone_reproduction_checks(pkg, checks)
    submission_process_docs_checks(pkg, checks)
    v04_manuscript_audit_checks(pkg, checks)
    reference_audit_checks(pkg, checks)
    manifest_checksum_checks(pkg, checks)

    write_report(report_path, zip_path, checks)
    print(report_path)
    failures = [check for check in checks if check.status != "PASS"]
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
