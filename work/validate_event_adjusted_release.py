#!/usr/bin/env python3
"""Validate the public event-adjusted manuscript release with the standard library."""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT = ARTICLE / "supplement"
CHINESE = ROOT / "outputs" / "cee_submission_latex_v0_8_chinese_review_article"
ADVISOR = ROOT / "outputs" / "advisor_chinese_manuscript_v1"
PERIODS = {0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_csv(name: str) -> list[dict[str, str]]:
    path = SUPPLEMENT / name
    require(path.is_file(), f"missing {path.relative_to(ROOT)}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def float_periods(rows: list[dict[str, str]]) -> set[float]:
    return {float(row["period_s"]) for row in rows}


def validate_files() -> None:
    required = [
        ARTICLE / "main.tex",
        ARTICLE / "main.pdf",
        ARTICLE / "supplementary_information.tex",
        ARTICLE / "supplementary_information.pdf",
        ARTICLE / "build_event_adjusted_supplement_figures.py",
        ROOT / "work" / "jshis_event_adjusted_station_model.py",
        ROOT / "work" / "audit_jshis_flatfile_selection.py",
        ROOT / "work" / "jshis_independent_gmpe_replication.py",
        ROOT / "work" / "jshis_station_uncertainty_propagation.py",
        ROOT / "work" / "jshis_spatial_model_complexity_audit.py",
        ROOT / "work" / "jshis_mf2013_applicability_audit.py",
        ROOT / "work" / "jshis_robustness_stress_tests.py",
        ROOT / "work" / "jshis_path_stratification_audit.py",
        ROOT / "work" / "jshis_balanced_station_model.py",
        ROOT / "work" / "verify_public_inputs.py",
        ROOT / "work" / "build_cee_initial_submission_package.py",
        ROOT / "work" / "validate_compact_peer_review_archive.py",
        ROOT / "environment.yml",
        ROOT / "public_inputs_manifest.tsv",
        ARTICLE / "figures" / "figure_ground_motion_model_sensitivity.pdf",
        ARTICLE / "figures" / "supplementary_figure_station_uncertainty.pdf",
        ARTICLE / "figures" / "supplementary_figure_robustness_stress_tests.pdf",
        ARTICLE / "figures" / "figure_path_stratification.pdf",
        ARTICLE / "figures" / "supplementary_figure_balanced_station_model.pdf",
        SUPPLEMENT / "jshis_flatfile_selection_audit.csv",
        SUPPLEMENT / "jshis_flatfile_selection_audit.md",
        SUPPLEMENT / "jshis_spatial_model_complexity_sensitivity.csv",
        SUPPLEMENT / "jshis_spatial_model_complexity_sensitivity.md",
        SUPPLEMENT / "jshis_mf2013_applicability_sensitivity.csv",
        SUPPLEMENT / "jshis_mf2013_applicability_sensitivity.md",
    ]
    for path in required:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty {path.relative_to(ROOT)}")

    with (ROOT / "public_inputs_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        inputs = list(csv.DictReader(handle, delimiter="\t"))
    require(
        {row["role"] for row in inputs}
        == {"strong_motion_flatfile", "mf2013_coefficients", "jshis_response_spectra"},
        "public input manifest roles are incomplete",
    )
    require(all(re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) for row in inputs), "invalid input SHA-256")
    require(all(row["source_url"].startswith("https://") for row in inputs), "input source URL is not HTTPS")


def validate_manuscript() -> None:
    text = (ARTICLE / "main.tex").read_text(encoding="utf-8")
    title_match = re.search(r"\\LARGE\\bfseries (.*?)\\par", text)
    require(title_match is not None, "title not found")
    title = title_match.group(1)
    title_words = re.findall(r"[A-Za-z0-9]+", title)
    require(len(title_words) <= 15, f"title has {len(title_words)} words")
    require(not re.search(r"[,.:;!?]", title), "title contains punctuation")

    abstract_source = text.split(r"\noindent\textbf{Abstract}", 1)[1].split(r"\noindent\textbf{Keywords", 1)[0]
    transition = re.search(r"Here we (?:analyse|present)", abstract_source)
    require(transition is not None, "abstract lacks a direct method transition")
    background = abstract_source[: transition.start()]
    require(background.count(".") == 2, "abstract does not contain two background sentences")
    abstract = re.sub(r"\\[A-Za-z]+(?:\{[^}]*\})?", " ", abstract_source)
    abstract_words = re.findall(r"[A-Za-z0-9]+(?:[.-][A-Za-z0-9]+)*", abstract)
    require(len(abstract_words) <= 150, f"abstract has {len(abstract_words)} words")
    require("222,664" in abstract and "surface records" in abstract, "abstract lacks the surface-record sample")
    for claim in ["RotD100", "0.890", "0.770", "12.6", "0.280", "1.338", "0.959", "9.9"]:
        require(claim in abstract, f"abstract lacks key result: {claim}")
    for claim in ["10.5\\%", "14.0\\%", "7.1\\%", "19.3\\%"]:
        require(claim in text, f"manuscript lacks model-complexity result: {claim}")

    subheadings = re.findall(r"\\subsection\{([^}]*)\}", text)
    require(all(len(heading) < 60 for heading in subheadings), "a Results or Methods subheading has 60 or more characters")
    require(
        all(not re.search(r"[,.:;!?]", heading) for heading in subheadings),
        "a Results or Methods subheading contains punctuation",
    )
    main_text = text.split(r"\section{Introduction}", 1)[1].split(r"\section{Methods}", 1)[0]
    main_text = re.sub(r"\\begin\{(?:figure|table)\}.*?\\end\{(?:figure|table)\}", " ", main_text, flags=re.S)
    main_text = re.sub(r"\\[A-Za-z]+(?:\[[^]]*\])?(?:\{[^}]*\})?", " ", main_text)
    main_words = re.findall(r"[A-Za-z0-9]+(?:[.-][A-Za-z0-9]+)*", main_text)
    require(len(main_words) <= 5_000, f"Introduction, Results and Discussion have {len(main_words)} words")
    display_items = len(re.findall(r"\\begin\{(?:figure|table)\}", text.split(r"\section{Methods}", 1)[0]))
    require(display_items <= 10, f"main manuscript has {display_items} display items")
    introduction = text.split(r"\section{Introduction}", 1)[1].split(r"\section{Results}", 1)[0]
    introduction_words = re.findall(r"[A-Za-z0-9]+(?:[.-][A-Za-z0-9]+)*", introduction)
    require(len(introduction_words) < 1_000, f"Introduction has {len(introduction_words)} words")
    last_introduction_paragraph = [paragraph for paragraph in introduction.split("\n\n") if paragraph.strip()][-1]
    for phrase in ["Event holdout", "path-stratified analyses", "average-path station adjustment"]:
        require(phrase in last_introduction_paragraph, f"Introduction scope paragraph lacks: {phrase}")
    require(
        "The manuscript and code were written by the authors. ChatGPT was used only for language and formatting revision and code verification" in text,
        "generative-AI assistance statement is not in the approved wording",
    )
    require(
        "Haoyu Zhou designed the analysis" in text and "Qiang Ma supervised the study" in text,
        "author contributions do not use full author names",
    )
    require("H.Z. designed" not in text and "Q.M. supervised" not in text, "initial-only author contributions remain")
    require("No generative-AI image is included" in text, "generative-image status is not disclosed")
    require(
        text.count("accompanying peer-review archive") >= 2,
        "data or code availability does not provide the peer-review archive",
    )
    require(
        "archived with a DOI before publication" in text,
        "code availability lacks the planned public DOI release",
    )
    require(
        "Analysis and figure-generation code is available at" not in text,
        "private GitHub repository is incorrectly described as currently public",
    )
    require(r"M_{\mathrm{JMA}}\geq5" in text, "JMA magnitude threshold is not stated")
    require("8,675 lack finite F-net $M_w$" in text, "finite-Mw sample loss is not stated")
    require(r"moment magnitude $M_w\geq5$" not in text, "public subset is incorrectly labelled as Mw>=5")
    require(
        "Each fold fits the two-way event-station decomposition separately" in text,
        "event holdout is not documented as separate two-way decompositions",
    )
    require("used only for manuscript-format checks" not in text, "obsolete AI-use statement remains")
    require(r"\section*{References}" not in text, "manual References heading duplicates the bibliography heading")
    for stale in ["322,020", "2,267", "23.1\\%", "0.589", "1.449", "0.048 g", "10,467"]:
        require(stale not in text, f"stale claim remains: {stale}")

    cover = (ARTICLE / "cover_letter_cee.md").read_text(encoding="utf-8")
    for claim in [
        "Repeatable station terms redistribute long-period response spectra across Japanese strong-motion sites",
        "Institute of Engineering Mechanics, China Earthquake Administration",
        "No. 29 Xuefu Road",
        "Harbin 150080",
        "maqiang@iem.ac.cn",
        "peer-review archive",
        "0.959",
        "9.9%",
        "0.649--1.445",
    ]:
        require(claim in cover, f"cover letter lacks: {claim}")


def validate_supplementary_order() -> None:
    main_text = (ARTICLE / "main.tex").read_text(encoding="utf-8")
    supplement_text = (ARTICLE / "supplementary_information.tex").read_text(encoding="utf-8")
    for kind, expected in [("Fig", 7), ("Table", 15)]:
        seen: list[int] = []
        for value in re.findall(rf"Supplementary {kind}\.?\s*(\d+)", main_text):
            number = int(value)
            if number not in seen:
                seen.append(number)
        require(seen == list(range(1, expected + 1)), f"Supplementary {kind} first-appearance order is {seen}")

    require(
        len(re.findall(r"\\begin\{table\}", supplement_text)) == 15,
        "Supplementary Information does not contain fifteen tables",
    )
    require(
        r"\renewcommand{\figurename}{Supplementary Figure}" in supplement_text
        and r"\renewcommand{\tablename}{Supplementary Table}" in supplement_text,
        "Supplementary figure or table naming is not configured",
    )
    require(r"\caption{Supplementary Table" not in supplement_text, "a supplementary table caption repeats its number")
    require(
        len(re.findall(r"\\begin\{figure\}", supplement_text)) == 7,
        "Supplementary Information does not contain seven figures",
    )
    require(
        r"\renewcommand{\refname}{Supplementary References}" in supplement_text
        and r"\section*{Supplementary References}" not in supplement_text,
        "Supplementary References heading is duplicated or mislabelled",
    )


def citation_order(text: str) -> tuple[list[str], list[str]]:
    body = text.split(r"\begin{thebibliography}", 1)[0]
    first_appearance: list[str] = []
    for group in re.findall(r"\\cite\{([^}]*)\}", body):
        for key in group.split(","):
            if key not in first_appearance:
                first_appearance.append(key)
    bibliography = re.findall(r"\\bibitem\{([^}]*)\}", text)
    return first_appearance, bibliography


def bibliography_entries(text: str) -> dict[str, str]:
    matches = list(re.finditer(r"\\bibitem\{([^}]*)\}", text))
    entries: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else text.find(r"\end{thebibliography}", match.end())
        entries[match.group(1)] = " ".join(text[match.end() : end].split())
    return entries


def validate_references() -> None:
    main_text = (ARTICLE / "main.tex").read_text(encoding="utf-8")
    supplement_text = (ARTICLE / "supplementary_information.tex").read_text(encoding="utf-8")
    for name, text, expected in [
        ("main manuscript", main_text, 33),
        ("Supplementary Information", supplement_text, 11),
    ]:
        cited, bibliography = citation_order(text)
        require(len(bibliography) == expected, f"unexpected {name} reference count: {len(bibliography)}")
        require(cited == bibliography, f"{name} references are missing, uncited or out of first-appearance order")

    shared = (ARTICLE / "references_shared.tex").read_text(encoding="utf-8")
    shared_keys = re.findall(r"\\bibitem\{([^}]*)\}", shared)
    main_keys = re.findall(r"\\bibitem\{([^}]*)\}", main_text)
    require(shared_keys == main_keys, "Chinese shared-reference keys differ from the English manuscript")
    require(
        bibliography_entries(shared) == bibliography_entries(main_text),
        "Chinese shared-reference content differs from the English manuscript",
    )
    main_entries = bibliography_entries(main_text)
    require(
        all("https://" in entry for entry in main_entries.values()),
        "a main reference lacks a DOI or official HTTPS source",
    )
    supplement_entries = bibliography_entries(supplement_text)
    require(
        all("https://" in entry for entry in supplement_entries.values()),
        "a supplementary reference lacks a DOI or official HTTPS source",
    )
    for name, text in [
        ("main manuscript", main_text),
        ("Chinese shared references", shared),
        ("Supplementary Information", supplement_text),
    ]:
        require("Dohi, Y. et al." in text, f"{name} does not abbreviate the six-author Dohi reference")
        require("Dohi, Y., Shigeno" not in text, f"{name} retains all six Dohi authors")


def validate_chinese_sync() -> None:
    if not (CHINESE / "main_zh.tex").is_file():
        return
    for path in [
        CHINESE / "main_zh.tex",
        CHINESE / "main_zh.pdf",
        CHINESE / "main_zh_preview_safe.pdf",
        ADVISOR / "main_zh.tex",
        ADVISOR / "main_zh.pdf",
        ADVISOR / "main_zh_preview_safe.pdf",
    ]:
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty local Chinese artifact: {path}")
    main_text = (ARTICLE / "main.tex").read_text(encoding="utf-8")
    chinese_text = (CHINESE / "main_zh.tex").read_text(encoding="utf-8")
    english_abstract = main_text.split(r"\noindent\textbf{Abstract}", 1)[1].split(
        r"\noindent\textbf{Keywords", 1
    )[0]
    chinese_english_abstract = chinese_text.split(r"\noindent{\bfseries Abstract:}", 1)[1].split(
        r"\noindent{\bfseries Keywords:", 1
    )[0]
    require(
        " ".join(english_abstract.split()) == " ".join(chinese_english_abstract.split()),
        "English abstract in the Chinese manuscript is not synchronized",
    )

    chinese_citations: list[str] = []
    for group in re.findall(r"\\cite\{([^}]*)\}", chinese_text):
        for key in group.split(","):
            if key not in chinese_citations:
                chinese_citations.append(key)
    main_bibliography = re.findall(r"\\bibitem\{([^}]*)\}", main_text)
    require(chinese_citations == main_bibliography, "Chinese citation order differs from the English bibliography")

    for stale in ["0.928", "0.780", "0.625", "1.366", "0.074 g", "1,153"]:
        require(stale not in chinese_text, f"stale Chinese claim remains: {stale}")
    for current in [
        "RotD100",
        "0.890",
        "0.770",
        "12.6\\%",
        "0.624",
        "1.338",
        "0.075 g",
        "0.959",
        "9.9\\%",
        "53,517",
        "488个地震",
        "$-2.6\\%$",
    ]:
        require(current in chinese_text, f"Chinese manuscript lacks current result: {current}")
    require("分别对4组训练事件和1组检验事件拟合完整的事件--台站双向固定效应模型" in chinese_text, "Chinese event-holdout method is stale")
    require(
        r"figure_path_stratification.pdf" in chinese_text,
        "Chinese manuscript lacks the path-stratification figure",
    )
    require("补充表12" in chinese_text and "补充图7" in chinese_text, "Chinese supplementary numbering is incomplete")
    require("补充表13" in chinese_text, "Chinese sample-flow and component table numbering is incomplete")
    require("补充表14" in chinese_text, "Chinese complexity-sensitivity table is not cited")
    require("补充表15" in chinese_text, "Chinese MF2013-applicability table is not cited")
    require("未复原Kanno-PGA距离截断" in chinese_text, "Chinese MF2013 applicability boundary is missing")
    require("整体中位数的变化方向随数据范围和台站集合而变" in chinese_text, "Chinese aggregate-median boundary is missing")
    require(r"M_{\mathrm{JMA}}\geq5" in chinese_text, "Chinese manuscript lacks the JMA magnitude threshold")
    require("8,675条缺少有限的F-net $M_w$" in chinese_text, "Chinese manuscript lacks the finite-Mw sample loss")
    require(chinese_text.count("随审稿档案提供") >= 2, "Chinese data or code availability is not synchronized")
    for claim in ["10.5\\%--14.0\\%", "7.1\\%--19.3\\%"]:
        require(claim in chinese_text, f"Chinese manuscript lacks model-complexity result: {claim}")


def validate_sample_selection() -> None:
    rows = read_csv("jshis_flatfile_selection_audit.csv")
    expected = [
        (1, "Public sub1-v2024 archive", 333_808, 0, 1_840, 2_581),
        (2, "Ground-surface installation", 231_380, 102_428, 1_840, 1_882),
        (3, "Finite F-net moment magnitude Mw", 222_705, 8_675, 1_737, 1_881),
        (4, "Supported source class and positive distance", 222_705, 0, 1_737, 1_881),
        (5, "Positive AVS30", 222_664, 41, 1_737, 1_880),
        (6, "Finite RotD100 and MF2013 predictions at all eight periods", 222_664, 0, 1_737, 1_880),
    ]
    observed = [
        (
            int(row["stage_order"]),
            row["stage"],
            int(row["records"]),
            int(row["excluded_at_stage"]),
            int(row["earthquakes"]),
            int(row["stations"]),
        )
        for row in rows
    ]
    require(observed == expected, f"flatfile sample flow differs from the audited release: {observed}")

    report = (SUPPLEMENT / "jshis_flatfile_selection_audit.md").read_text(encoding="utf-8")
    for claim in [
        "JMA magnitude range: 5.0--9.0",
        "Shortest-fault-distance range: 1.0001--299.9953 km",
        "Surface records with nonpositive AVS30: 42",
        "also lack finite Mw: 1",
        "3 s: 222,664",
    ]:
        require(claim in report, f"flatfile selection report lacks: {claim}")


def validate_mf2013_applicability() -> None:
    rows = read_csv("jshis_mf2013_applicability_sensitivity.csv")
    require(len(rows) == 8 and float_periods(rows) == PERIODS, "MF2013 applicability period coverage mismatch")
    require(all(float(row["minimum_mw"]) == 5.5 for row in rows), "MF2013 sensitivity magnitude screen changed")
    require(
        all(float(row["maximum_distance_km_exclusive"]) == 200.0 for row in rows),
        "MF2013 sensitivity distance screen changed",
    )
    require(
        all(int(float(row["minimum_event_stations"])) == 5 for row in rows),
        "MF2013 sensitivity event-station screen changed",
    )
    require(all(int(float(row["n_records"])) == 53_517 for row in rows), "MF2013 sensitivity record count changed")
    require(all(int(float(row["n_events"])) == 488 for row in rows), "MF2013 sensitivity event count changed")
    require(
        all(int(float(row["n_eligible_stations"])) == 772 for row in rows),
        "MF2013 sensitivity eligible-station count changed",
    )
    require(
        max(float(row["solver_max_parameter_change"]) for row in rows) <= 1.01e-10,
        "MF2013 sensitivity decomposition did not converge",
    )

    sa3 = next(row for row in rows if float(row["period_s"]) == 3.0)
    bounds = {
        "station_field_correlation": (0.9593, 0.9595),
        "station_field_difference_q95_log10": (0.0827, 0.0830),
        "restricted_event_holdout_correlation": (0.8946, 0.8948),
        "restricted_spatial_rmse_gain_pct": (9.8, 9.9),
        "oof_prediction_correlation": (0.7702, 0.7704),
        "restricted_multiplier_q05": (0.648, 0.650),
        "restricted_multiplier_q50": (0.963, 0.965),
        "restricted_multiplier_q95": (1.444, 1.446),
        "common_surface_sa_median_g": (0.0890, 0.0892),
        "common_primary_adjusted_sa_median_g": (0.0910, 0.0912),
        "common_restricted_adjusted_sa_median_g": (0.0876, 0.0879),
    }
    for column, (low, high) in bounds.items():
        value = float(sa3[column])
        require(low < value < high, f"unexpected MF2013 applicability {column}: {value}")
    require(int(float(sa3["common_hazard_stations"])) == 772, "MF2013 common hazard-station count changed")

    short = next(row for row in rows if float(row["period_s"]) == 0.1)
    require(float(short["restricted_spatial_rmse_gain_pct"]) < 0, "adverse 0.1 s applicability result disappeared")
    require(
        all(float(row["restricted_spatial_rmse_gain_pct"]) > 0 for row in rows if float(row["period_s"]) >= 0.2),
        "a 0.2--5.0 s applicability spatial result is not positive",
    )


def validate_model_complexity() -> None:
    rows = read_csv("jshis_spatial_model_complexity_sensitivity.csv")
    require(len(rows) == 288, f"unexpected model-complexity metric count: {len(rows)}")
    require({int(row["spatial_seed"]) for row in rows} == {20_260_710}, "complexity audit seed changed")
    settings = {
        "primary",
        "fewer_leaves",
        "more_leaves",
        "larger_minimum_leaf",
        "stronger_l2",
        "shorter_iteration_budget",
    }
    require({row["setting"] for row in rows} == settings, "model-complexity settings are incomplete")
    overall = [row for row in rows if row["scope"] == "overall"]
    require(len(overall) == 48 and float_periods(overall) == PERIODS, "complexity overall coverage mismatch")
    gains = [float(row["rmse_reduction_vs_zero_pct"]) for row in overall]
    require(7.0 < min(gains) < 7.2 and 19.2 < max(gains) < 19.4, "complexity gain range changed")
    sa3 = [row for row in overall if float(row["period_s"]) == 3.0]
    sa3_gains = [float(row["rmse_reduction_vs_zero_pct"]) for row in sa3]
    require(10.4 < min(sa3_gains) < 10.5 and 13.9 < max(sa3_gains) < 14.1, "SA3 complexity range changed")

    primary = {
        (float(row["period_s"]), row["scope"], int(row["fold"])): row
        for row in rows
        if row["setting"] == "primary"
    }
    released = {
        (float(row["period_s"]), row["scope"], int(row["fold"])): row
        for row in read_csv("jshis_event_adjusted_station_model_metrics.csv")
        if row["model"] == "physical_spatial_hgb"
    }
    require(primary.keys() == released.keys(), "complexity primary rows differ from released spatial metrics")
    for key in primary:
        for field in [
            "weighted_rmse_log10",
            "zero_baseline_rmse_log10",
            "rmse_reduction_vs_zero_pct",
            "observed_predicted_correlation",
            "prediction_centering_shift_log10",
        ]:
            require(
                abs(float(primary[key][field]) - float(released[key][field])) < 1e-9,
                f"complexity primary does not reproduce released metric {key} {field}",
            )


def validate_decomposition() -> None:
    rows = read_csv("jshis_event_station_decomposition_metrics.csv")
    unique = {(float(row["period_s"]), row["model"]) for row in rows}
    require(len(unique) == 24, f"expected 24 decompositions, found {len(unique)}")
    require(float_periods(rows) == PERIODS, "decomposition period coverage mismatch")
    maximum_change = max(float(row["max_parameter_change"]) for row in rows)
    require(maximum_change <= 1.01e-10, f"fixed effects not converged: {maximum_change}")

    terms = read_csv("jshis_event_adjusted_station_terms.csv")
    require(len(terms) == 45_120, f"unexpected station-term row count: {len(terms)}")
    require(
        {int(float(row["installation_situation_id"])) for row in terms} == {1},
        "non-surface station terms entered the primary analysis",
    )
    require({row["response_component"] for row in terms} == {"RotD100"}, "primary station terms are not RotD100")
    weighted_sum: dict[tuple[float, str], float] = defaultdict(float)
    weight: dict[tuple[float, str], float] = defaultdict(float)
    for row in terms:
        key = (float(row["period_s"]), row["model"])
        count = float(row["n_records"])
        value = float(row["station_effect_log10"])
        require(math.isfinite(value), "non-finite station term")
        weighted_sum[key] += count * value
        weight[key] += count
    largest_mean = max(abs(weighted_sum[key] / weight[key]) for key in weighted_sum)
    require(largest_mean < 1e-12, f"station terms are not zero-centred: {largest_mean}")


def validate_prediction() -> None:
    rows = read_csv("jshis_event_adjusted_station_model_metrics.csv")
    overall = [row for row in rows if row["scope"] == "overall" and row["model"] == "physical_spatial_hgb"]
    require(float_periods(overall) == PERIODS, "model period coverage mismatch")
    sa3 = next(row for row in overall if float(row["period_s"]) == 3.0)
    gain = float(sa3["rmse_reduction_vs_zero_pct"])
    correlation = float(sa3["observed_predicted_correlation"])
    require(12.5 < gain < 12.7, f"unexpected SA3 RMSE gain: {gain}")
    require(0.49 < correlation < 0.51, f"unexpected SA3 correlation: {correlation}")
    folds = [
        row
        for row in rows
        if row["scope"] == "fold" and row["model"] == "physical_spatial_hgb" and float(row["period_s"]) == 3.0
    ]
    require(len(folds) == 5, "missing SA3 spatial folds")
    require(min(float(row["rmse_reduction_vs_zero_pct"]) for row in folds) > 0, "a SA3 spatial fold does not improve")

    repeat = read_csv("jshis_event_holdout_station_repeatability.csv")
    require(len(repeat) == 48, f"unexpected event-holdout row count: {len(repeat)}")
    require({row["decomposition_method"] for row in repeat} == {"two_way_fixed_effects"}, "event holdout is not two-way")
    require({int(float(row["event_split_seed"])) for row in repeat} == {20_260_710}, "event split seed changed")
    require(
        max(float(row["train_max_parameter_change"]) for row in repeat) <= 1.01e-10
        and max(float(row["test_max_parameter_change"]) for row in repeat) <= 1.01e-10,
        "an event-holdout fixed-effect fit did not converge",
    )
    repeat3 = next(row for row in repeat if row["scope"] == "fold_mean" and float(row["period_s"]) == 3.0)
    require(0.88 < float(repeat3["train_test_station_correlation"]) < 0.90, "unexpected event-holdout correlation")
    require(67.0 < float(repeat3["rmse_reduction_vs_zero_pct"]) < 68.5, "unexpected event-holdout gain")


def validate_hazard() -> None:
    values = read_csv("jshis_event_adjusted_surface_spectrum_values.csv")
    require(len(values) == 52_096, f"unexpected response-spectrum row count: {len(values)}")
    value_keys = {
        (row["siteid2"], float(row["period_s"]), row["probability_level"])
        for row in values
    }
    require(len(value_keys) == len(values), "duplicate station-period-probability response-spectrum rows")

    rows = read_csv("jshis_event_adjusted_surface_spectrum_summary.csv")
    require(float_periods(rows) == PERIODS, "response-spectrum period coverage mismatch")
    sa3 = next(
        row
        for row in rows
        if float(row["period_s"]) == 3.0 and row["probability_level"] == "50y_10pct"
    )
    bounds = {
        "oof_multiplier_q05": (0.62, 0.63),
        "oof_multiplier_q50": (0.88, 0.90),
        "oof_multiplier_q95": (1.33, 1.35),
        "ergodic_surface_sa_g_q50": (0.080, 0.082),
        "adjusted_surface_sa_g_q50": (0.073, 0.076),
    }
    for column, (low, high) in bounds.items():
        value = float(sa3[column])
        require(low < value < high, f"unexpected {column}: {value}")


def validate_independent_gmpe() -> None:
    terms = read_csv("jshis_zhao2006_station_terms.csv")
    require(len(terms) == 15_040, f"unexpected Zhao station-term row count: {len(terms)}")
    require(
        {int(float(row["installation_situation_id"])) for row in terms} == {1},
        "non-surface records entered the Zhao station terms",
    )
    require({row["response_component"] for row in terms} == {"RotD50"}, "Zhao station terms are not RotD50")
    replication = read_csv("jshis_independent_gmpe_station_term_replication.csv")
    require(len(replication) == 8 and float_periods(replication) == PERIODS, "Zhao replication period coverage mismatch")
    sa3 = next(row for row in replication if float(row["period_s"]) == 3.0)
    require(0.77 < float(sa3["pearson_correlation"]) < 0.79, "unexpected MF2013-Zhao SA3 correlation")
    require(int(float(sa3["n_paired_stations"])) == 1_628, "unexpected Zhao paired-station count")

    repeat = read_csv("jshis_zhao2006_event_holdout_repeatability.csv")
    require(len(repeat) == 48, f"unexpected Zhao event-holdout row count: {len(repeat)}")
    require({row["decomposition_method"] for row in repeat} == {"two_way_fixed_effects"}, "Zhao event holdout is not two-way")
    require({int(float(row["event_split_seed"])) for row in repeat} == {20_260_710}, "Zhao event split seed changed")
    require(
        max(float(row["train_max_parameter_change"]) for row in repeat) <= 1.01e-10
        and max(float(row["test_max_parameter_change"]) for row in repeat) <= 1.01e-10,
        "a Zhao event-holdout fixed-effect fit did not converge",
    )
    repeat3 = next(row for row in repeat if row["scope"] == "fold_mean" and float(row["period_s"]) == 3.0)
    require(0.94 < float(repeat3["train_test_station_correlation"]) < 0.96, "unexpected Zhao event repeatability")

    primary_repeat = read_csv("jshis_event_holdout_station_repeatability.csv")
    primary_folds = {
        (float(row["period_s"]), int(float(row["fold"]))): (
            int(float(row["event_split_seed"])),
            int(float(row["n_train_events"])),
            int(float(row["n_test_events"])),
            int(float(row["n_paired_stations"])),
            int(float(row["test_record_weight"])),
        )
        for row in primary_repeat
        if row["scope"] == "fold"
    }
    zhao_folds = {
        (float(row["period_s"]), int(float(row["fold"]))): (
            int(float(row["event_split_seed"])),
            int(float(row["n_train_events"])),
            int(float(row["n_test_events"])),
            int(float(row["n_paired_stations"])),
            int(float(row["test_record_weight"])),
        )
        for row in repeat
        if row["scope"] == "fold"
    }
    require(primary_folds == zhao_folds, "MF2013 and Zhao event-holdout folds differ")
    metrics = read_csv("jshis_zhao2006_station_model_metrics.csv")
    model3 = next(
        row
        for row in metrics
        if row["scope"] == "overall"
        and row["model"] == "physical_spatial_hgb"
        and float(row["period_s"]) == 3.0
    )
    require(37.0 < float(model3["rmse_reduction_vs_zero_pct"]) < 38.5, "unexpected Zhao SA3 spatial gain")
    require(0.78 < float(model3["observed_predicted_correlation"]) < 0.81, "unexpected Zhao SA3 spatial correlation")


def validate_response_component() -> None:
    rows = read_csv("jshis_response_component_audit.csv")
    require(len(rows) == 8 and float_periods(rows) == PERIODS, "response-component period coverage mismatch")
    require(
        all(int(float(row["n_records_rotd100"])) == 222_664 for row in rows),
        "RotD100 response-component sample changed",
    )
    require(
        all(int(float(row["n_records_rotd50"])) == 222_664 for row in rows),
        "RotD50 response-component sample changed",
    )
    sa3 = next(row for row in rows if float(row["period_s"]) == 3.0)
    require(0.9985 < float(sa3["station_pearson"]) < 0.9987, "unexpected SA3 component correlation")
    require(
        0.0070 < float(sa3["station_difference_weighted_rmse_log10"]) < 0.0072,
        "unexpected SA3 component difference",
    )
    require(
        0.817 < float(sa3["record_rotd50_over_rotd100_q50"]) < 0.821,
        "unexpected SA3 RotD50/RotD100 ratio",
    )


def validate_path_stratification() -> None:
    rows = read_csv("jshis_path_stratification_stability.csv")
    require(len(rows) == 88 and float_periods(rows) == PERIODS, "path-stratification coverage mismatch")
    require(
        all(int(float(row["solver_nonconverged_components"])) == 0 for row in rows),
        "a path-stratum component did not converge",
    )
    require(
        max(float(row["solver_max_relative_normal_residual"]) for row in rows) < 1e-8,
        "path-stratum normal-equation residual is too large",
    )
    require(
        all("weighted_within_component_correlation" in row for row in rows),
        "path-stratum component-invariant correlation is missing",
    )

    sa3 = [row for row in rows if float(row["period_s"]) == 3.0]
    azimuth = [row for row in sa3 if row["stratification"] == "hypocentral_azimuth"]
    source = [row for row in sa3 if row["stratification"] == "source_type"]
    require(len(azimuth) == 4 and len(source) == 3, "SA3 path-stratum rows are incomplete")
    require(
        min(float(row["weighted_station_term_correlation"]) for row in azimuth) < 0.30,
        "low directional full-field agreement disappeared",
    )
    require(
        min(float(row["weighted_station_term_correlation"]) for row in source) > 0.92,
        "source-class station fields are no longer stable",
    )

    split = read_csv("jshis_path_stratification_split_half_repeatability.csv")
    require(len(split) == 165, f"unexpected path split-half row count: {len(split)}")
    require(float_periods(split) == {1.0, 2.0, 3.0}, "path split-half period coverage mismatch")
    require(
        all("weighted_split_half_correlation_within_component" in row for row in split),
        "component-invariant split-half correlation is missing",
    )
    directional = [
        row
        for row in split
        if float(row["period_s"]) == 3.0
        and row["stratification"] == "hypocentral_azimuth"
        and row["stratum"] == "000-090 deg"
    ]
    require(len(directional) == 5, "directional split-half repeats are incomplete")
    mean_directional = sum(
        float(row["weighted_split_half_correlation_within_component"]) for row in directional
    ) / len(directional)
    require(0.88 < mean_directional < 0.89, "directional within-component repeatability changed")
    sa3_azimuth = [
        row
        for row in split
        if float(row["period_s"]) == 3.0 and row["stratification"] == "hypocentral_azimuth"
    ]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in sa3_azimuth:
        grouped[row["stratum"]].append(row)
    within_means = {
        stratum: sum(float(row["weighted_split_half_correlation_within_component"]) for row in block)
        / len(block)
        for stratum, block in grouped.items()
    }
    require(0.763 < min(within_means.values()) < 0.765, "unexpected minimum SA3 sector repeatability")
    require(0.968 < max(within_means.values()) < 0.970, "unexpected maximum SA3 sector repeatability")
    require(
        max(
            abs(
                float(row["weighted_split_half_correlation"])
                - float(row["weighted_split_half_correlation_within_component"])
            )
            for row in sa3_azimuth
        )
        < 0.04,
        "SA3 directional repeatability is sensitive to component alignment",
    )
    require(
        max(int(row["n_component_pairs"]) for row in split) > 1,
        "split-half component-pair audit did not exercise disconnected cases",
    )

    metrics = read_csv("jshis_balanced_station_model_metrics.csv")
    require(len(metrics) == 288 and float_periods(metrics) == PERIODS, "equal-stratum metric coverage mismatch")
    selected = [
        row
        for row in metrics
        if row["model"] == "physical_spatial_hgb"
        and row["scope"] == "overall"
        and float(row["period_s"]) in {1.0, 2.0, 3.0}
    ]
    require(len(selected) == 9, "equal-stratum 1--3 s summary is incomplete")
    source3 = next(
        row
        for row in selected
        if row["residual_model"] == "source_balanced" and float(row["period_s"]) == 3.0
    )
    require(8.9 < float(source3["rmse_reduction_vs_zero_pct"]) < 9.1, "unexpected source-balanced SA3 gain")
    balanced_folds = [
        row
        for row in metrics
        if row["model"] == "physical_spatial_hgb"
        and row["scope"] == "fold"
        and float(row["period_s"]) in {1.0, 2.0, 3.0}
    ]
    require(
        min(float(row["rmse_reduction_vs_zero_pct"]) for row in balanced_folds) < -60,
        "adverse equal-stratum spatial fold disappeared",
    )


def validate_uncertainty() -> None:
    intervals = read_csv("jshis_station_term_prediction_intervals.csv")
    require(len(intervals) == 13_024, f"unexpected station interval row count: {len(intervals)}")
    keys = {(row["siteid2"], float(row["period_s"])) for row in intervals}
    require(len(keys) == len(intervals), "duplicate station-period prediction intervals")
    require(
        all(float(row["lower_prediction_log10"]) <= float(row["point_prediction_log10"]) <= float(row["upper_prediction_log10"]) for row in intervals),
        "a station prediction interval does not contain its point estimate",
    )

    coverage = read_csv("jshis_station_term_interval_coverage.csv")
    overall = [row for row in coverage if row["scope"] == "overall"]
    require(float_periods(overall) == PERIODS, "uncertainty period coverage mismatch")
    sa3 = next(row for row in overall if float(row["period_s"]) == 3.0)
    require(88.5 < float(sa3["station_coverage_pct"]) < 89.5, "unexpected SA3 station interval coverage")
    require(88.5 < float(sa3["record_weighted_coverage_pct"]) < 89.5, "unexpected weighted SA3 coverage")
    require(2.8 < float(sa3["median_interval_factor_span"]) < 3.1, "unexpected SA3 interval span")

    hazard = read_csv("jshis_surface_spectrum_prediction_intervals.csv")
    require(len(hazard) == 52_096, f"unexpected propagated interval row count: {len(hazard)}")
    hazard_keys = {(row["siteid2"], float(row["period_s"]), row["probability_level"]) for row in hazard}
    require(len(hazard_keys) == len(hazard), "duplicate propagated spectrum intervals")
    require(
        all(
            float(row["adjusted_surface_sa_lower_g"])
            <= float(row["adjusted_surface_sa_point_g"])
            <= float(row["adjusted_surface_sa_upper_g"])
            for row in hazard
        ),
        "a propagated spectrum interval does not contain its point estimate",
    )


def validate_stress_tests() -> None:
    metrics = read_csv("jshis_station_transfer_stress_metrics.csv")
    require(len(metrics) == 168, f"unexpected transfer metric row count: {len(metrics)}")
    common_network = [
        row
        for row in metrics
        if row["validation"] == "network_transfer" and row["model"] == "common_feature_spatial_hgb"
    ]
    require(len(common_network) == 16 and float_periods(common_network) == PERIODS, "common-feature network coverage mismatch")
    require(min(float(row["rmse_reduction_vs_zero_pct"]) for row in common_network) > 0, "common-feature network transfer is not positive in every direction and period")

    full_region3 = [
        row
        for row in metrics
        if row["validation"] == "macroregion_leaveout"
        and row["model"] == "physical_spatial_hgb"
        and float(row["period_s"]) == 3.0
    ]
    require(len(full_region3) == 5, "missing SA3 macroregion tests")
    require(min(float(row["rmse_reduction_vs_zero_pct"]) for row in full_region3) > 0, "full SA3 model fails a macroregion")
    common_region3 = [
        row
        for row in metrics
        if row["validation"] == "macroregion_leaveout"
        and row["model"] == "common_feature_spatial_hgb"
        and float(row["period_s"]) == 3.0
    ]
    require(min(float(row["rmse_reduction_vs_zero_pct"]) for row in common_region3) < 0, "adverse common-feature regional result disappeared")

    event = read_csv("jshis_influential_event_stability.csv")
    require(len(event) == 18, f"unexpected event-stability row count: {len(event)}")
    single = [row for row in event if row["perturbation"] == "single_event_removal"]
    joint = [row for row in event if row["perturbation"].startswith("joint_top10")]
    require(len(single) == 10 and len(joint) == 8, "influential-event test coverage mismatch")
    require(min(float(row["weighted_station_term_correlation"]) for row in single) > 0.9999, "single-event stability degraded")
    require(min(float(row["weighted_station_term_correlation"]) for row in joint) > 0.9989, "joint-event stability degraded")


def main() -> None:
    validate_files()
    validate_manuscript()
    validate_supplementary_order()
    validate_references()
    validate_chinese_sync()
    validate_sample_selection()
    validate_mf2013_applicability()
    validate_decomposition()
    validate_prediction()
    validate_model_complexity()
    validate_hazard()
    validate_independent_gmpe()
    validate_response_component()
    validate_path_stratification()
    validate_uncertainty()
    validate_stress_tests()
    print("event-adjusted release validation passed")


if __name__ == "__main__":
    main()
