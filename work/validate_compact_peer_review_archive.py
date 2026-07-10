#!/usr/bin/env python3
"""Validate the compact CEE peer-review archive using only the standard library."""

from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "outputs" / "cee_submission_latex_v0_8_english_article"
SUPPLEMENT = ARTICLE / "supplement"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(name: str) -> list[dict[str, str]]:
    with (SUPPLEMENT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def verify_manifest() -> None:
    manifest = ROOT / "MANIFEST_SHA256.tsv"
    require(manifest.is_file(), "MANIFEST_SHA256.tsv is missing")
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(len(rows) >= 50, "review archive manifest is unexpectedly short")
    for row in rows:
        path = ROOT / row["path"]
        require(path.is_file(), f"manifest file is missing: {row['path']}")
        require(path.stat().st_size == int(row["bytes"]), f"size mismatch: {row['path']}")
        require(sha256(path) == row["sha256"], f"SHA-256 mismatch: {row['path']}")


def validate_claim_tables() -> None:
    selection = read_csv("jshis_flatfile_selection_audit.csv")
    require([int(row["records"]) for row in selection] == [333_808, 231_380, 222_705, 222_705, 222_664, 222_664], "record-selection flow changed")

    decomposition = read_csv("jshis_event_station_decomposition_metrics.csv")
    unique_fits = {(float(row["period_s"]), row["model"]) for row in decomposition}
    require(len(unique_fits) == 24, "decomposition coverage changed")
    require(max(float(row["max_parameter_change"]) for row in decomposition) < 1e-10, "a decomposition did not converge")

    event = [row for row in read_csv("jshis_event_holdout_station_repeatability.csv") if row["scope"] == "fold_mean"]
    event3 = next(row for row in event if float(row["period_s"]) == 3.0)
    require(0.889 < float(event3["train_test_station_correlation"]) < 0.891, "SA3 event repeatability changed")
    require(67.6 < float(event3["rmse_reduction_vs_zero_pct"]) < 67.8, "SA3 event RMSE gain changed")

    spatial = [
        row
        for row in read_csv("jshis_event_adjusted_station_model_metrics.csv")
        if row["model"] == "physical_spatial_hgb" and row["scope"] == "overall"
    ]
    spatial3 = next(row for row in spatial if float(row["period_s"]) == 3.0)
    require(12.5 < float(spatial3["rmse_reduction_vs_zero_pct"]) < 12.7, "SA3 spatial gain changed")
    require(0.499 < float(spatial3["observed_predicted_correlation"]) < 0.501, "SA3 spatial correlation changed")

    zhao3 = next(row for row in read_csv("jshis_independent_gmpe_station_term_replication.csv") if float(row["period_s"]) == 3.0)
    require(0.769 < float(zhao3["pearson_correlation"]) < 0.771, "SA3 MF2013-Zhao agreement changed")

    path = [row for row in read_csv("jshis_path_stratification_stability.csv") if float(row["period_s"]) == 3.0 and row["stratification"] == "hypocentral_azimuth"]
    path_values = [float(row["weighted_station_term_correlation"]) for row in path]
    require(0.279 < min(path_values) < 0.281 and 0.885 < max(path_values) < 0.887, "SA3 directional range changed")

    split = [row for row in read_csv("jshis_path_stratification_split_half_repeatability.csv") if float(row["period_s"]) == 3.0 and row["stratification"] == "hypocentral_azimuth"]
    groups: dict[str, list[float]] = defaultdict(list)
    for row in split:
        groups[row["stratum"]].append(float(row["weighted_split_half_correlation_within_component"]))
    means = [sum(values) / len(values) for values in groups.values()]
    require(0.763 < min(means) < 0.765 and 0.968 < max(means) < 0.970, "SA3 within-component repeatability changed")

    complexity = [row for row in read_csv("jshis_spatial_model_complexity_sensitivity.csv") if row["scope"] == "overall"]
    require(len(complexity) == 48 and min(float(row["rmse_reduction_vs_zero_pct"]) for row in complexity) > 7.0, "model-complexity sensitivity changed")

    applicability = read_csv("jshis_mf2013_applicability_sensitivity.csv")
    require(len(applicability) == 8, "MF2013 applicability coverage changed")
    require(all(float(row["minimum_kanno_pga_cm_s2"]) == 10.0 for row in applicability), "Kanno-PGA cutoff changed")
    require(all(int(float(row["n_records"])) == 35_857 for row in applicability), "MF2013 applicability record count changed")
    require(all(int(float(row["n_events"])) == 411 for row in applicability), "MF2013 applicability event count changed")
    applicability3 = next(row for row in applicability if float(row["period_s"]) == 3.0)
    require(0.9593 < float(applicability3["station_field_correlation"]) < 0.9595, "SA3 applicability field correlation changed")
    require(4.1 < float(applicability3["restricted_spatial_rmse_gain_pct"]) < 4.2, "SA3 applicability spatial gain changed")
    require(0.659 < float(applicability3["restricted_multiplier_q05"]) < 0.661, "SA3 restricted lower multiplier changed")
    require(1.502 < float(applicability3["restricted_multiplier_q95"]) < 1.504, "SA3 restricted upper multiplier changed")
    applicability01 = next(row for row in applicability if float(row["period_s"]) == 0.1)
    require(float(applicability01["restricted_spatial_rmse_gain_pct"]) < 0, "adverse 0.1 s applicability result disappeared")

    kiknet = read_csv("jshis_kiknet_surface_borehole_validation.csv")
    require(len(kiknet) == 56, "KiK-net paired-sensor coverage changed")
    kiknet3 = next(
        row for row in kiknet if float(row["period_s"]) == 3.0 and row["scope"] == "full"
    )
    require(int(float(kiknet3["n_paired_records"])) == 102_428, "KiK-net paired count changed")
    require(0.314 < float(kiknet3["station_transfer_pearson"]) < 0.316, "SA3 paired-sensor correlation changed")
    kiknet3_mean = next(
        row for row in kiknet if float(row["period_s"]) == 3.0 and row["scope"] == "fold_mean"
    )
    require(0.289 < float(kiknet3_mean["test_station_term_pearson"]) < 0.291, "SA3 held-event paired prediction changed")
    require(2.7 < float(kiknet3_mean["transfer_prediction_rmse_gain_pct"]) < 2.8, "SA3 held-event paired gain changed")

    hazard3 = next(
        row
        for row in read_csv("jshis_event_adjusted_surface_spectrum_summary.csv")
        if float(row["period_s"]) == 3.0 and row["probability_level"] == "50y_10pct"
    )
    for field, lower, upper in [
        ("official_vs400_sa_g_q50", 0.077, 0.079),
        ("ergodic_surface_sa_g_q50", 0.080, 0.082),
        ("adjusted_surface_sa_g_q50", 0.074, 0.076),
        ("oof_multiplier_q05", 0.623, 0.625),
        ("oof_multiplier_q95", 1.337, 1.339),
    ]:
        require(lower < float(hazard3[field]) < upper, f"SA3 hazard summary changed: {field}")

    coverage3 = next(
        row
        for row in read_csv("jshis_station_term_interval_coverage.csv")
        if float(row["period_s"]) == 3.0 and row["scope"] == "overall"
    )
    require(88.9 < float(coverage3["station_coverage_pct"]) < 89.1, "SA3 interval coverage changed")
    require(2.9 < float(coverage3["median_interval_factor_span"]) < 3.1, "SA3 interval span changed")

    component3 = next(row for row in read_csv("jshis_response_component_audit.csv") if float(row["period_s"]) == 3.0)
    require(0.9985 < float(component3["station_pearson"]) < 0.9987, "SA3 component correlation changed")
    require(0.0179 < float(component3["station_difference_abs_q95_log10"]) < 0.0181, "SA3 component difference changed")


def validate_manuscript_sources() -> None:
    main = (ARTICLE / "main.tex").read_text(encoding="utf-8")
    supplement = (ARTICLE / "supplementary_information.tex").read_text(encoding="utf-8")
    for claim in [
        "222,664",
        "0.890",
        "0.770",
        "102,428",
        "0.315",
        "0.290",
        "12.6\\%",
        "0.280--0.886",
        "0.624--1.338",
        "0.959",
        "4.1\\%",
        "35,857",
        "$-2.3\\%$",
        "0.660--1.503",
    ]:
        require(claim in main, f"main manuscript lacks: {claim}")
    require(main.count("accompanying peer-review archive") >= 2, "review archive is not declared")
    require(supplement.count(r"\begin{table}") == 16, "Supplementary Information table count changed")
    require(supplement.count(r"\begin{figure}") == 7, "Supplementary Information figure count changed")


def main() -> None:
    verify_manifest()
    validate_claim_tables()
    validate_manuscript_sources()
    print("compact peer-review archive validation passed")


if __name__ == "__main__":
    main()
