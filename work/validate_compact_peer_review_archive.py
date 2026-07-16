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

    shape = read_csv("jshis_kiknet_spectral_shape_validation.csv")
    require(len(shape) == 7, "KiK-net spectral-shape coverage changed")
    shape_full = next(row for row in shape if row["scope"] == "full")
    shape_mean = next(row for row in shape if row["scope"] == "held_event_fold_mean")
    require(0.676 < float(shape_full["shape_pearson"]) < 0.678, "full spectral-shape correlation changed")
    require(0.635 < float(shape_mean["shape_pearson"]) < 0.636, "held-event spectral-shape correlation changed")

    transfer = [
        row
        for row in read_csv("jshis_cross_network_transfer_metrics.csv")
        if row["scope"] == "station_aggregate"
        and row["source_network"] == "K-NET"
        and row["target_network"] == "KiK-net"
        and row["model"] == "physical_spatial_hgb"
    ]
    require(len(transfer) == 8, "cross-network period coverage changed")
    transfer3 = next(row for row in transfer if float(row["period_s"]) == 3.0)
    require(0.597 < float(transfer3["pearson"]) < 0.599, "SA3 cross-network correlation changed")
    require(20.7 < float(transfer3["rmse_gain_pct"]) < 20.8, "SA3 cross-network gain changed")

    temporal = [
        row
        for row in read_csv("jshis_temporal_network_transfer_metrics.csv")
        if float(row["train_fraction"]) == 0.8
        and row["source_network"] == "K-NET"
        and row["target_network"] == "KiK-net"
        and row["model"] == "physical_spatial_hgb"
    ]
    require(len(temporal) == 8, "chronological transfer period coverage changed")
    require(
        min(float(row["pearson_ci_low"]) for row in temporal) > 0
        and min(float(row["rmse_gain_pct_ci_low"]) for row in temporal) > 0,
        "a chronological transfer interval crosses zero",
    )
    temporal3 = next(row for row in temporal if float(row["period_s"]) == 3.0)
    require(0.535 < float(temporal3["pearson"]) < 0.537, "SA3 chronological correlation changed")
    require(19.3 < float(temporal3["rmse_gain_pct"]) < 19.5, "SA3 chronological gain changed")
    reverse_temporal3 = next(
        row
        for row in read_csv("jshis_temporal_network_transfer_metrics.csv")
        if float(row["period_s"]) == 3.0
        and float(row["train_fraction"]) == 0.8
        and row["source_network"] == "KiK-net"
        and row["target_network"] == "K-NET"
        and row["model"] == "physical_spatial_hgb"
    )
    require(float(reverse_temporal3["rmse_gain_pct"]) < 0, "adverse chronological reverse result disappeared")

    temporal_splits = read_csv("jshis_temporal_network_event_splits.csv")
    require(len(temporal_splits) == 30_648, "chronological event-split coverage changed")
    split_groups: dict[tuple[float, float], list[dict[str, str]]] = defaultdict(list)
    for row in temporal_splits:
        split_groups[(float(row["period_s"]), float(row["train_fraction"]))].append(row)
    require(len(split_groups) == 24, "chronological event-split groups changed")
    for rows in split_groups.values():
        early = [row["origin_time"] for row in rows if row["partition"] == "early_train"]
        late = [row["origin_time"] for row in rows if row["partition"] == "late_test"]
        require(early and late and max(early) < min(late), "chronological event split leaks across time")

    esm_selection = read_csv("esm_external_selection_audit.csv")
    esm_primary_selection = next(
        row
        for row in esm_selection
        if row["analysis"] == "bindi2014_rhyp_mw4.5-7.6_2015-01-01"
        and row["stage"] == "iterative_event_station_support"
    )
    require(int(esm_primary_selection["n_records"]) == 13_430, "ESM primary record count changed")
    require(int(esm_primary_selection["n_events"]) == 634, "ESM primary event count changed")
    esm_means = [
        row
        for row in read_csv("esm_external_repeatability_metrics.csv")
        if row["scope"] == "fold_mean"
        and row["backbone"] == "bindi2014_rhyp"
        and float(row["minimum_mw"]) == 4.5
    ]
    require(len(esm_means) == 3, "ESM primary period coverage changed")
    esm3 = next(row for row in esm_means if float(row["period_s"]) == 3.0)
    require(0.9277 < float(esm3["pearson"]) < 0.9280, "ESM SA3 recurrence changed")
    require(66.3 < float(esm3["rmse_gain_pct"]) < 66.5, "ESM SA3 gain changed")
    require(float(esm3["pearson_ci_low"]) > 0.89, "ESM SA3 recurrence interval degraded")
    esm_transfer3 = next(
        row
        for row in read_csv("esm_cross_region_transfer_metrics.csv")
        if float(row["period_s"]) == 3.0
    )
    require(
        esm_transfer3["prediction_centering"] == "target_prediction_weighted_zero_without_labels",
        "Japan-to-ESM centring rule changed",
    )
    require(0.046 < float(esm_transfer3["pearson"]) < 0.049, "Japan-to-ESM SA3 transfer changed")
    require(float(esm_transfer3["rmse_gain_pct"]) < 0, "adverse Japan-to-ESM transfer disappeared")

    esm_spatial = read_csv("esm_local_spatial_prediction_metrics.csv")
    require(
        {int(float(row["n_spatial_blocks"])) for row in esm_spatial} == {3},
        "ESM spatial-block count changed",
    )
    esm_spatial3 = next(
        row
        for row in esm_spatial
        if row["scope"] == "overall"
        and row["model"] == "common_site_hgb"
        and float(row["period_s"]) == 3.0
    )
    esm_event_space3 = next(
        row
        for row in esm_spatial
        if row["scope"] == "event_fold_mean"
        and row["model"] == "common_site_hgb"
        and float(row["period_s"]) == 3.0
    )
    require(0.073 < float(esm_spatial3["pearson"]) < 0.075, "ESM local SA3 correlation changed")
    require(-3.8 < float(esm_spatial3["rmse_gain_pct"]) < -3.5, "ESM local SA3 gain changed")
    require(0.020 < float(esm_event_space3["pearson"]) < 0.022, "ESM event-space SA3 correlation changed")
    require(-3.9 < float(esm_event_space3["rmse_gain_pct"]) < -3.5, "ESM event-space SA3 gain changed")
    require(
        all(
            float(row["rmse_gain_pct"]) < 0
            for row in esm_spatial
            if row["scope"] in {"overall", "event_fold_mean"}
        ),
        "an adverse ESM spatial result disappeared",
    )

    source_summary = read_csv("jshis_source_category_surface_spectrum_summary.csv")
    require(len(source_summary) == 96, "source-category summary coverage changed")
    source3 = {
        row["source_category"]: row
        for row in source_summary
        if float(row["period_s"]) == 3.0 and row["probability_level"] == "50y_10pct"
    }
    require(set(source3) == {"all_earthquakes", "active_shallow", "subduction"}, "source-category SA3 rows missing")
    require(0.030 < float(source3["active_shallow"]["official_vs400_sa_g_q50"]) < 0.032, "active-shallow median changed")
    require(0.070 < float(source3["subduction"]["official_vs400_sa_g_q50"]) < 0.072, "subduction median changed")
    source_stations = read_csv("jshis_source_category_sa3_station_comparison.csv")
    source_counts: dict[str, int] = defaultdict(int)
    for row in source_stations:
        source_counts[row["dominant_source_category"]] += 1
    require(
        dict(source_counts) == {"subduction": 1_428, "active_shallow": 200},
        f"source-category station counts changed: {dict(source_counts)}",
    )

    conditioned = [
        row
        for row in read_csv("jshis_source_conditioned_station_metrics.csv")
        if row["model"] == "physical_spatial_hgb" and row["scope"] == "overall"
    ]
    require(len(conditioned) == 24, "source-conditioned metric coverage changed")
    conditioned_by_key = {
        (row["source_type"], float(row["period_s"])): row for row in conditioned
    }
    require(
        14.2
        < float(conditioned_by_key[("Intraplate", 3.0)]["rmse_reduction_vs_zero_pct"])
        < 14.4
        and 14.7
        < float(conditioned_by_key[("Intraplate", 5.0)]["rmse_reduction_vs_zero_pct"])
        < 14.9,
        "source-conditioned intraplate result changed",
    )
    require(
        float(conditioned_by_key[("Crustal", 0.1)]["rmse_reduction_vs_zero_pct"]) < 0
        and float(conditioned_by_key[("Interplate", 2.0)]["rmse_reduction_vs_zero_pct"]) < 0,
        "an adverse source-conditioned result disappeared",
    )

    sung = read_csv("sung2025_kanto_sa5_external_metrics.csv")
    require(len(sung) == 6, "Sung comparison coverage changed")
    sung_primary = next(row for row in sung if row["comparison_field"] == "crustal_observed")
    require(int(sung_primary["n_stations"]) == 364, "Sung matched count changed")
    require(0.303 < float(sung_primary["pearson"]) < 0.305, "Sung Pearson result changed")
    require(0.16 < float(sung_primary["pearson_ci_low"]) < 0.17, "Sung lower interval changed")

    interval = read_csv("jshis_hazard_impact_interval_robustness.csv")
    interval3 = next(
        row
        for row in interval
        if float(row["period_s"]) == 3.0 and float(row["confidence_level_pct"]) == 90.0
    )
    require(1.9 < float(interval3["robust_direction_fraction_pct"]) < 2.1, "SA3 interval robustness changed")
    agreement3 = next(
        row
        for row in read_csv("jshis_hazard_impact_independent_model_summary.csv")
        if float(row["period_s"]) == 3.0
    )
    require(0.735 < float(agreement3["independent_model_pearson"]) < 0.737, "SA3 independent-field correlation changed")
    require(79.6 < float(agreement3["direction_agreement_pct"]) < 79.7, "SA3 direction agreement changed")

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
        "0.598",
        "20.8\\%",
        "0.536",
        "19.4\\%",
        "1,277",
        "0.635",
        "0.304",
        "0.163--0.434",
        "14.3\\%",
        "14.8\\%",
        "79.7\\%",
        "2.0\\%",
        "12.6\\%",
        "0.280--0.886",
        "0.624--1.338",
        "0.959",
        "4.1\\%",
        "35,857",
        "$-2.3\\%$",
        "0.660--1.503",
        "13,430",
        "0.928",
        "0.047",
        "$-3.6\\%$",
        "$-7.6\\%$",
        "1,428 of 1,628",
    ]:
        require(claim in main, f"main manuscript lacks: {claim}")
    require(main.count("accompanying peer-review archive") >= 2, "review archive is not declared")
    require(supplement.count(r"\begin{table}") == 26, "Supplementary Information table count changed")
    require(supplement.count(r"\begin{figure}") == 14, "Supplementary Information figure count changed")


def main() -> None:
    verify_manifest()
    validate_claim_tables()
    validate_manuscript_sources()
    print("compact peer-review archive validation passed")


if __name__ == "__main__":
    main()
