# CEE manuscript package v0.8

Date: 2026-06-18.

This directory contains the single-column Communications Earth & Environment submission package. The text uses the verified J-SHIS/NIED, MF2013, residual-analysis, held-out station prediction, official response-spectrum sensitivity, public source-parameter checks, and independent site-variable audit results.

The PDF is an A4 single-column manuscript with 11 pt body text and 1.20 line spacing. Figures use the English publication figures from the project output directory.

The manuscript now includes author affiliations, correspondence, ORCID identifiers, author contributions, no-funding statement, competing-interests statement, and GitHub code/data availability. References follow numbered Nature style for CEE: in-text superscript numbers and an ordered reference list.

Figure 1 is a six-panel study overview and audit summary. The main text now uses nine figures, including spatial evidence, split official-spectrum sensitivity figures, representative spectra, and uncertainty checks. The Supplementary Information contains Supplementary Figs. S1--S15 and Supplementary Tables S1--S27.

The 2026-06-18 submission-closure revision tightens the abstract, opening argument, discussion, data/code wording, cover letter, official-source boundary, and claim-evidence traceability table. The paper is framed around one contribution: Japanese strong-motion records contain a 1--3 s non-ergodic station term that is predictable from public site variables and large enough to change official long-period response-spectrum coordinates.

The latest refresh adds public source-parameter checks and an interpretable SA(3.0 s) site-variable model. The manuscript cites GitHub as the current release location and leaves Zenodo DOI minting as a final pre-submission or journal-requested step.

The latest text pass clarifies four submission-facing points in the main manuscript: the concrete contribution relative to residual studies, the role and boundary of the response-spectrum sensitivity calculation, the path from held-out station terms to spatially continuous site terms, and the external geologic information needed for stronger physical attribution.

The OpenQuake/GEM Japan-model PSHA sensitivity check has been added as Supplementary Fig. S4 and Supplementary Table S8. It is framed as an independent engine-level sensitivity check. The 2026-06-15 update adds the 5,322-mesh public-parameter audit, five-point sigma sensitivity audit, residual-root-cause audit, 10,467-mesh expanded execution audit, 10,467-mesh near-station physical audit, basin-period proxy audit, tail-robustness case-profile audit, nearest-station official borehole-profile audit, sampled-grid SA(3.0 s) continuous correction-surface audit, independent GSI terrain geomorphic proxy audit, and J-SHIS VS400 site-amplification audit as main-text evidence and Supplementary Figs. S5--S12 / Supplementary Tables S9--S19.

The main manuscript uses nine main figures. The Supplementary Information documents the external-model check, model-extension diagnostics, OpenQuake PSHA sensitivity check, sigma sensitivity audit, residual-tail audit, basin-period proxy audit, tail-robustness case-profile audit, nearest-station borehole-profile audit, sampled-grid correction-surface audit, GSI terrain geomorphic proxy audit, J-SHIS VS400 site-amplification audit, multi-event KiK-net transfer-function audit, all-period response-spectrum propagation audit, official PSHA input-chain audit, public source-parameter checks, and interpretable SA(3.0 s) site-variable model.

## Files

- `main.tex`: LaTeX source file.
- `main.pdf`: compiled English manuscript.
- `supplementary_information.tex`: LaTeX source for Supplementary Figs. S1--S15 and Supplementary Tables S1--S27.
- `supplementary_information.pdf`: compiled Supplementary Information file.
- `cover_letter_cee.md`: concise Communications Earth & Environment cover letter.
- `reference_full_metadata_audit_2026-06-14.md`: full reference metadata audit against DOI, publisher, and J-SHIS/NIED sources.
- `../cee_v08_public_archive_readiness.md`: pre-deposition public-archive readiness review for the v0.8 package.
- `../cee_v08_public_archive_file_review.csv`: file-level archive review table.
- `../cee_release_metadata_zenodo_v0_8_draft.json`: Zenodo metadata file to finalize after author approval.
- `figures/figure1.pdf`--`figures/figure9_engineering_cases_uncertainty.pdf`: English main figure PDFs.
- `figures/figure9_engineering_cases_uncertainty.pdf` / `.png`: main-text representative city-nearest spectra and uncertainty figure.
- `figures/figure5_spatial_evidence.pdf` / `.png`: main-text spatial evidence figure linking sampled-grid correction, gridded AVS30, validation, and response-spectrum sensitivity.
- `figures/openquake_jpn_station_psha_sensitivity.pdf` / `.png`: Supplementary Fig. S4, independent OpenQuake/GEM Japan-model PSHA sensitivity check.
- `figures/figure10_sigma_sensitivity.png`: Supplementary Fig. S5, five-point sigma sensitivity audit.
- `figures/figure11_residual_root_cause.png`: Supplementary Fig. S6, 5,322-mesh residual-tail diagnostics.
- `figures/figure12_basin_period_proxy.png`: Supplementary Fig. S7, 10,467-mesh basin-period proxy audit.
- `figures/figure13_tail_robustness_cases.png`: Supplementary Fig. S8, 10,467-mesh tail-robustness and case-profile audit.
- `figures/figure14_nearest_station_borehole_audit.png`: Supplementary Fig. S9, 10,467-mesh nearest-station official borehole-profile audit.
- `figures/figure15_sa3_continuous_surface.png`: Supplementary Fig. S10, sampled-grid SA(3.0 s) continuous correction-surface audit.
- `figures/figure16_gsi_terrain_geomorphic_audit.png`: Supplementary Fig. S11, independent GSI terrain geomorphic proxy audit.
- `figures/figure17_zamp_vs400_site_amplification_audit.png`: Supplementary Fig. S12, J-SHIS VS400 gridded site-amplification audit.
- `build_figure1_en.py`: reproducible script for the six-panel overview and audit figure.
- `build_figure5_spatial_evidence_en.py`: reproducible script for the main-text spatial evidence figure.
- `build_figure8_en.py`: reproducible script for the official-data SA(3.0 s) response-spectrum check.
- `build_figure9_en.py`: reproducible script for representative city-nearest spectra and uncertainty checks.
- `supplement/jshis_nonergodic_station_model.md`: station-held-out residual correction report.
- `supplement/jshis_nonergodic_station_model_summary.csv`: summary metrics by target, split, and model.
- `supplement/jshis_nonergodic_station_model_fold_metrics.csv`: fold-level validation metrics.
- `supplement/jshis_station_model_extensions.md`: station-model extension report.
- `supplement/jshis_station_model_ablation_summary.csv`: spatial-block feature-set ablation summary.
- `supplement/jshis_station_model_ablation_fold_metrics.csv`: fold-level ablation metrics.
- `supplement/jshis_station_model_grouped_permutation_importance.csv`: grouped permutation importance metrics.
- `supplement/jshis_station_model_network_transfer.csv`: K-NET/KiK-net transfer metrics.
- `supplement/jshis_station_engineering_examples.csv`: SA(3.0 s) engineering-scale station correction examples.
- `supplement/jshis_station_hazard_impact.md`: SA(3.0 s) station-multiplier report used before the official-data sensitivity check.
- `supplement/jshis_station_hazard_impact_by_station.csv`: station-level SA(3.0 s) multipliers.
- `supplement/jshis_station_hazard_impact_summary.csv`: station-multiplier distribution summary.
- `supplement/jshis_station_hazard_impact_examples.csv`: example stations for engineering-scale interpretation.
- `supplement/jshis_station_hazard_curve_points.csv`: legacy curve-source table retained for audit traceability.
- `supplement/jshis_official_hazard_response_check.md`: official J-SHIS hazard-curve API and response-spectrum map check.
- `supplement/jshis_official_hazard_response_summary.csv`: official-data sensitivity-check summary.
- `supplement/jshis_official_response_sa3_station_values.csv`: matched station-level official and station-corrected SA(3.0 s) values.
- `supplement/jshis_official_response_uhs_examples.csv`: official and SA(3.0 s)-adjusted UHS examples.
- `supplement/jshis_official_hazard_curve_api_examples.csv`: API-derived example hazard curves.
- `supplement/jshis_official_hazard_response_fold_summary.csv`: spatial-block fold summary for the official-data sensitivity check.
- `supplement/jshis_official_hazard_response_bootstrap_ci.csv`: bootstrap intervals for official-data sensitivity-check quantiles.
- `supplement/jshis_official_hazard_response_uncertainty_spatial_check.md`: spatial and uncertainty check report.
- `supplement/jshis_station_region_leaveout_summary.csv`: seven-region geographic leave-out validation summary.
- `supplement/jshis_station_region_leaveout_fold_metrics.csv`: region-level leave-out metrics.
- `supplement/jshis_station_variable_substitution_summary.csv`: official-variable substitution summary.
- `supplement/jshis_station_variable_substitution_fold_metrics.csv`: fold-level variable-substitution metrics.
- `supplement/jshis_station_site_variable_collinearity_pairs.csv`: Spearman collinearity checks among public site variables.
- `supplement/jshis_station_site_variable_vif.csv`: variance-inflation factors for VS20, AVS30, D1400, and Dbase.
- `supplement/jshis_positive_sa3_station_correction_summary.csv`: positive SA(3.0 s) correction count and distribution summary.
- `supplement/jshis_positive_sa3_station_correction_audit.csv`: station-level positive-correction audit table.
- `supplement/jshis_station_matching_missingness_audit.csv`: record/site matching and missingness audit.
- `supplement/jshis_official_data_provenance_2026-06-14.csv`: official data-source provenance table.
- `supplement/cee_code_data_release_manifest_2026-06-14.md`: release manifest for derived data and code.
- `supplement/cee_claim_evidence_traceability.md` and `supplement/cee_claim_evidence_traceability.csv`: claim-evidence traceability files mapping the main claims to figures, derived tables, numerical anchors, interpretation boundaries, and release status.
- `supplement/cee_release_metadata_zenodo_draft.json`: Zenodo metadata file with DOI fields to complete after DOI minting.
- `supplement/jshis_existing_sensitivity_inventory.csv`: inventory of event, station, magnitude, distance, and external-model sensitivity outputs.
- `supplement/jshis_existing_influence_and_external_model_sensitivity.md`: Markdown summary of existing sensitivity outputs.
- `supplement/openquake_jpn_station_psha_summary.csv`: representative-station SA(3.0 s) PSHA sensitivity summary.
- `supplement/openquake_jpn_station_psha_uhs.csv`: baseline and station-shifted OpenQuake UHS ordinates.
- `supplement/openquake_jpn_station_psha_hazard_curves.csv`: baseline and station-shifted SA(3.0 s) hazard-curve source data.
- `supplement/openquake_jpn_station_psha_sensitivity.md`: OpenQuake/GEM Japan-model PSHA sensitivity report.
- `supplement/jshis_public_psha_national_stratified5322_arrays12_recompute_audit.md`: 5,322-mesh public-parameter national stratified audit.
- `supplement/jshis_public_psha_national_stratified5322_arrays12_run_summary.csv`: shard and tile timing summary for the 5,322-mesh audit.
- `supplement/jshis_public_psha_sigma_sensitivity1797_audit.md`: five-point sigma sensitivity audit.
- `supplement/jshis_public_psha_sigma_sensitivity1797_metrics.csv`: sigma sensitivity residual metrics.
- `supplement/jshis_public_psha_sigma_sensitivity1797_best.csv`: probability-level best sigma summary.
- `supplement/jshis_public_psha_residual_root_cause5322_audit.md`: residual-root-cause audit for the 5,322-mesh run.
- `supplement/jshis_public_psha_residual_root_cause5322_by_probability.csv`: residual metrics by official probability level.
- `supplement/jshis_public_psha_residual_root_cause5322_by_coarse_geo.csv`: residual metrics by coarse geographic band.
- `supplement/jshis_public_psha_residual_root_cause5322_by_official_bv_quantile.csv`: residual metrics by official PBV decile.
- `supplement/jshis_public_psha_residual_root_cause5322_by_candidate_quantile.csv`: residual metrics by candidate-source-count decile.
- `supplement/jshis_public_psha_residual_root_cause5322_spearman.csv`: monotonic associations with residual and absolute residual.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_audit.md`: 10,467-mesh public-parameter national stratified execution audit.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_run_summary.csv`: shard and tile timing summary for the 10,467-mesh audit.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_by_probability.csv`: residual metrics by official PBV probability level for the 10,467-mesh audit.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_by_coarse_geo.csv`: residual metrics by coarse geographic band for the 10,467-mesh audit.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_by_official_bv_quantile.csv`: residual metrics by official PBV decile for the 10,467-mesh audit.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_by_candidate_quantile.csv`: residual metrics by candidate-source-count decile for the 10,467-mesh audit.
- `supplement/jshis_public_psha_national_stratified10467_arrays8_recompute_spearman.csv`: monotonic associations with residual and absolute residual for the 10,467-mesh audit.
- `supplement/jshis_public_psha_10467_near_station_physical_audit.md`: near-station physical audit summary for the 10,467-mesh public-parameter run.
- `supplement/jshis_public_psha_10467_near_station_physical_summary.csv`: group-level basin-proxy and station-residual contrasts for the near-station physical audit.
- `supplement/jshis_public_psha_10467_near_station_physical_link.csv`: mesh-level nearest-station link table for the 10,467-mesh physical audit.
- `supplement/jshis_public_psha_10467_basin_period_proxy_audit.md`: basin-period proxy audit summary for the 10,467-mesh public-parameter run.
- `supplement/jshis_public_psha_10467_basin_period_proxy_summary.csv`: group-level D1400-derived period proxy and residual-tail metrics.
- `supplement/jshis_public_psha_10467_basin_period_proxy_associations.csv`: monotonic association and rank-test metrics for the basin-period proxy audit.
- `supplement/jshis_public_psha_10467_basin_period_proxy_tail_risk.csv`: residual-tail risk ratios for basin-period and geographic contrasts.
- `supplement/jshis_public_psha_10467_basin_period_proxy_link.csv`: mesh-level basin-period proxy link table.
- `supplement/jshis_public_psha_10467_tail_robustness_cases_audit.md`: tail-robustness and matched-case audit summary.
- `supplement/jshis_public_psha_10467_tail_robustness_distance_bins.csv`: distance-bin residual-tail summaries.
- `supplement/jshis_public_psha_10467_tail_robustness_stratified.csv`: stratified residual-tail risk ratios.
- `supplement/jshis_public_psha_10467_tail_threshold_sensitivity.csv`: residual-tail threshold sensitivity for 0.50--0.90 absolute log10 units.
- `supplement/jshis_public_psha_10467_tail_robustness_logistic.csv`: adjusted logistic tail-model coefficients.
- `supplement/jshis_public_psha_10467_tail_case_profiles.csv`: matched tail and non-tail mesh profiles.
- `supplement/jshis_public_psha_10467_tail_case_probability_rows.csv`: probability-level residual rows for matched case profiles.
- `supplement/jshis_public_psha_10467_nearest_station_borehole_audit.md`: nearest-station official borehole-profile audit summary.
- `supplement/jshis_public_psha_10467_nearest_station_borehole_summary.csv`: group-level borehole VS20, low-velocity thickness, and impedance-depth proxy contrasts.
- `supplement/jshis_public_psha_10467_nearest_station_borehole_profiles.csv`: station-level parsed public NIED borehole-profile metrics.
- `supplement/jshis_public_psha_10467_nearest_station_borehole_mesh_link.csv`: mesh-level link table joining the 10,467 sampled meshes to nearest-station borehole metrics.
- `supplement/jshis_sa3_continuous_correction_surface_audit.md`: sampled-grid SA(3.0 s) continuous correction-surface audit summary.
- `supplement/jshis_sa3_continuous_correction_surface_summary.csv`: validation and mesh-group summary for the correction-surface audit.
- `supplement/jshis_sa3_continuous_correction_surface_station_validation.csv`: spatial-block held-out station predictions for the correction-surface audit.
- `supplement/jshis_sa3_continuous_correction_surface_mesh.csv`: 10,467-mesh sampled-grid correction-surface values.
- `supplement/jshis_gsi_terrain_geomorphic_audit.md`: independent GSI terrain geomorphic proxy audit summary.
- `supplement/jshis_gsi_terrain_geomorphic_summary.csv`: group-level terrain summary for the GSI audit.
- `supplement/jshis_gsi_terrain_geomorphic_associations.csv`: association and rank-test metrics for the GSI audit.
- `supplement/jshis_gsi_terrain_geomorphic_mesh.csv`: mesh-level GSI elevation, slope, relief, and lowland-flat fields.
- `supplement/jshis_zamp_vs400_site_amplification_audit.md`: J-SHIS VS400 site-amplification audit summary.
- `supplement/jshis_zamp_vs400_site_amplification_summary.csv`: group-level AVS30 and amplification summary for the VS400 audit.
- `supplement/jshis_zamp_vs400_site_amplification_associations.csv`: association, rank-test, and tail-risk metrics for the VS400 audit.
- `supplement/jshis_zamp_vs400_site_amplification_jcode.csv`: engineering-geomorphology code summaries for the VS400 audit.
- `supplement/jshis_zamp_vs400_site_amplification_mesh.csv`: mesh-level VS400 site-amplification fields joined to the 10,467-mesh audit grid.
- `supplement/kiknet_surface_downhole_audit.md`: earlier paired KiK-net surface/downhole waveform proxy audit summary retained for comparison.
- `supplement/kiknet_surface_downhole_ratios.csv`: earlier record-level paired KiK-net surface/downhole amplitude and spectral-ratio proxy table retained for comparison.
- `supplement/kiknet_surface_downhole_audit.py`: script used to recompute the earlier paired KiK-net surface/downhole proxy audit from matched miniSEED records.
- `supplement/kiknet_multievent_transfer_function_audit.md`: two-event KiK-net surface/downhole transfer-function audit summary.
- `supplement/kiknet_multievent_transfer_functions.csv`: event-station transfer-function rows for 460 KiK-net surface/downhole pairs.
- `supplement/kiknet_multievent_transfer_function_bins.csv`: frequency-bin transfer ratios for the two-event KiK-net audit.
- `supplement/kiknet_multievent_transfer_function_event_summary.csv`: event-level transfer-function summary.
- `supplement/kiknet_multievent_transfer_function_station_summary.csv`: station-level transfer-function summary.
- `supplement/kiknet_multievent_transfer_function_audit.py`: script used to recompute the two-event KiK-net transfer-function audit from local KiK-net waveform inputs.
- `supplement/jshis_multiperiod_nonergodic_psha_audit.md`: all-period official response-spectrum station-correction audit summary.
- `supplement/jshis_multiperiod_nonergodic_psha_summary.csv`: station and sampled-grid response-spectrum correction summaries.
- `supplement/jshis_multiperiod_station_corrections.csv`: station correction anchors at SA(0.3 s), SA(1.0 s), and SA(3.0 s).
- `supplement/jshis_multiperiod_station_surface_validation.csv`: spatial-block validation rows for the all-period correction surface.
- `supplement/jshis_multiperiod_station_surface_validation_summary.csv`: validation summary for the all-period correction surface.
- `supplement/jshis_multiperiod_continuous_correction_surface_mesh.csv`: 10,467-mesh correction surface at response-spectrum periods.
- `supplement/jshis_multiperiod_official_response_station_values.csv`: official and corrected response-spectrum ordinates at matched station meshes.
- `supplement/jshis_multiperiod_official_response_mesh_values.csv`: official and corrected response-spectrum ordinates at the sampled audit meshes.
- `supplement/jshis_multiperiod_nonergodic_psha_audit.py`: script used to recompute the all-period official-map propagation audit from local J-SHIS response-spectrum inputs.
- `supplement/jshis_official_psha_input_chain_audit.md`: official J-SHIS 2024 PSHM input-chain audit summary.
- `supplement/jshis_official_psha_input_audit_summary.csv`: category-level summary for official PSHA input files.
- `supplement/jshis_official_psha_input_file_inventory.csv`: CSV parameter-file inventory for the official input-chain audit.
- `supplement/jshis_official_psha_input_shape_inventory.csv`: shapefile inventory for the official input-chain audit.
- `supplement/jshis_official_psha_input_download_manifest.csv`: download manifest for official PSHA parameter packages.
- `supplement/jshis_official_psha_input_audit.py`: script used to recompute the official input-chain audit from local J-SHIS parameter-package inputs.
- `supplement/jshis_representative_city_spectrum_cases.csv`: representative city-nearest official and station-corrected response spectra for Fig. 9.
- `supplement/jshis_representative_uncertainty_summary.csv`: bootstrap and spatial-block uncertainty summary for Fig. 9.
- `render/contact_sheet.png`: rendered-page overview for local visual layout checking; excluded from the clean zip package.

## Added experiment

- `jshis_nonergodic_station_residual_model.py` in the project `work/` directory runs the station-held-out residual correction experiment.
- `jshis_nonergodic_station_model_extensions.py` in the project `work/` directory runs the ablation, grouped permutation importance, network-transfer, engineering-example, and reproducibility-pack steps.
- `jshis_station_hazard_impact.py` in the project `work/` directory builds the station-multiplier audit tables.
- `jshis_official_hazard_response_check.py` in the project `work/` directory applies the station multipliers to official J-SHIS 2020 response-spectrum map ordinates at matched 250 m station meshes and stores API-derived reference curves.
- `jshis_official_hazard_response_uncertainty_map.py` in the project `work/` directory builds spatial-block, bootstrap, and map-based checks for Figure 9.
- `jshis_cee_additional_experiments.py` in the project `work/` directory builds regional leave-out, official-variable substitution, positive-correction, matching/missingness, provenance, and release-manifest audit outputs.
- `build_openquake_jpn_station_psha.py` and `postprocess_openquake_jpn_station_psha.py` in the project `work/` directory build and post-process the independent OpenQuake/GEM Japan-model PSHA sensitivity check.
- `jshis_public_psha_tiled_map_recompute.py`, `jshis_public_psha_national_stratified_sample.py`, and related `jshis_public_psha_*.py` scripts in the project `work/` directory build the public-parameter J-SHIS PBV audit, national stratified samples, tiled recomputation outputs, sigma sensitivity checks, residual-tail diagnostics, near-station physical audit, basin-period proxy audit, tail-robustness case-profile audit, and nearest-station official borehole-profile audit.
- `jshis_official_psha_input_audit.py` in the project `work/` directory audits the official J-SHIS 2024 PSHM input-chain packages and defines the source-model provenance boundary.
- `jshis_sa3_continuous_correction_surface.py` in the project `work/` directory builds the sampled-grid SA(3.0 s) continuous correction-surface audit.
- `jshis_multiperiod_nonergodic_psha_audit.py` in the project `work/` directory propagates station-correction terms through official J-SHIS response-spectrum ordinates at 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, and 5.0 s.
- `build_figure9_en.py` in this article directory builds representative city-nearest spectra and uncertainty summaries from the derived response-spectrum tables.
- `kiknet_multievent_transfer_function_audit.py` in the project `work/` directory builds the two-event KiK-net surface/downhole transfer-function audit.
- `jshis_gsi_terrain_geomorphic_audit.py` in the project `work/` directory builds the independent GSI terrain geomorphic proxy audit.
- `jshis_zamp_vs400_site_amplification_audit.py` in the project `work/` directory builds the J-SHIS VS400 gridded site-amplification audit.
- `audit_cee_v08_public_archive.py` in the project `work/` directory builds the v0.8 public-archive readiness review and Zenodo metadata file.
- The experiment uses MF2013 site-model station residuals with at least 20 records per station and target.
- It reports random station holdout and spatial-block holdout reductions in weighted station-residual RMSE.
- The extension reports spatial-block ablations, feature-group importance, K-NET/KiK-net transfer, and SA(3.0 s) correction factors.
- The official-response-spectrum sensitivity check propagates station-correction terms through the available official J-SHIS response-spectrum ordinates. The official input-chain audit documents the public source-model packages, and the calculation keeps source recurrence, path terms, and logic-tree aggregation in the official product.
- The uncertainty check resamples matched stations 2000 times and reports spatial-block fold summaries.
- The OpenQuake sensitivity check uses packaged OpenQuake/GEM Japan mosaic inputs for five representative stations and applies the deterministic SA(3.0 s) station multiplier to the computed hazard curves and UHS ordinates.
- The public-parameter J-SHIS PBV audit reports a 10,467-mesh national stratified execution audit, 10,467-mesh residual-stratification tables, a 10,467-mesh near-station physical audit, a basin-period proxy audit, a tail-robustness case-profile audit, a nearest-station official borehole-profile audit, an official PSHA input-chain audit, sampled-grid correction-surface audits across SA(0.3 s), SA(1.0 s), and SA(3.0 s), an all-period official response-spectrum propagation audit, an independent GSI terrain geomorphic proxy audit, and a J-SHIS VS400 gridded site-amplification audit in the manuscript package. The 5,322-mesh residual-tail files remain for detailed comparison with the earlier diagnostic subset.

## Compile

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

## Remaining items before submission

- Confirm that the GitHub repository is public: `https://github.com/zhouhaoyiu/japan-strong-motion-site-terms`.
- Create a Zenodo DOI from the GitHub release if the journal or coauthors require a DOI-bearing archive in addition to GitHub.
- Refresh J-SHIS/NIED web access dates on the actual submission date.
- Re-export the final bibliography through the target journal system or a reference manager; Morikawa et al. (2024) should be checked against the final WCEE proceedings record if page or paper-number metadata become available.
- Apply the final target journal template after the journal is fixed.
