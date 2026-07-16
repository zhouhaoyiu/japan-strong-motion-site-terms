# CEE Manuscript Package

## Manuscript

- `main.tex` and `main.pdf`: CEE Article manuscript.
- `supplementary_information.tex` and `supplementary_information.pdf`: Supplementary Information.
- `cover_letter_cee.md`: cover-letter draft.
- `advisor_summary_cn.md`: Chinese project summary for advisor review.
- `reference_traceability.md`: reference metadata and local full-text audit.

## Current figures

The main manuscript uses eight figures:

1. `figure_event_adjusted_overview`
2. `figure_event_adjusted_site_structure`
3. `figure_ground_motion_model_sensitivity`
4. `figure_cross_network_transfer_validation`
5. `figure_esm_external_validation`
6. `figure_event_adjusted_multiperiod`
7. `figure_path_stratification`
8. `figure_event_adjusted_city_cases`

The Supplementary Information uses fourteen figures: event repeatability, national paired-sensor amplitude and spectral-shape validation, spatial validation, robustness stress tests, equal-stratum prediction, source-conditioned fields, the independent Kanto 3D-simulation comparison, response-spectrum sensitivity, station-model uncertainty, independent correction-field stability, official source-category response sensitivity, the KiK-net waveform comparison and chronological cross-network transfer. Vector PDFs are the manuscript sources; PNG copies are retained for portal compatibility.

## Derived tables

The `supplement/` directory contains the primary MF2013 results, the RotD50--RotD100 component audit, Zhao 2006 model sensitivity, the complete MF2013 regression-domain screen, national KiK-net paired-sensor validation, random-event and chronological network transfer, the post-2015 ESM recurrence and spatial tests, station-model prediction intervals, independent correction-field stability, source-conditioned station fields, the Kanto 3D-simulation comparison, official J-SHIS source-category response sensitivity, influential-event deletion results, path/source stratification and equal-stratum sensitivity. Rebuild the analyses in this order:

```bash
conda run -n japan-station-terms python ../../work/verify_public_inputs.py
conda run -n japan-station-terms python ../../work/audit_jshis_flatfile_selection.py
conda run -n japan-station-terms python ../../work/jshis_event_adjusted_station_model.py
conda run -n japan-station-terms python ../../work/jshis_spatial_model_complexity_audit.py
conda run -n japan-station-terms python ../../work/jshis_mf2013_applicability_audit.py
conda run -n japan-station-terms python ../../work/jshis_kiknet_surface_borehole_validation.py
conda run -n japan-station-terms python ../../work/jshis_cross_network_transfer_validation.py
conda run -n japan-station-terms python ../../work/jshis_temporal_network_transfer_validation.py
conda run -n japan-station-terms python ../../work/esm_external_station_validation.py
conda run -n japan-station-terms python ../../work/jshis_source_category_response_sensitivity.py
conda run -n japan-station-terms python ../../work/jshis_source_conditioned_station_fields.py
conda run -n japan-station-terms python ../../work/sung2025_kanto_external_comparison.py
conda run -n japan-station-terms python ../../work/jshis_independent_gmpe_replication.py
conda run -n japan-station-terms python ../../work/jshis_station_uncertainty_propagation.py
conda run -n japan-station-terms python ../../work/jshis_hazard_impact_robustness.py
conda run -n japan-station-terms python ../../work/jshis_robustness_stress_tests.py
conda run -n japan-station-terms python ../../work/jshis_path_stratification_audit.py
conda run -n japan-station-terms python ../../work/jshis_balanced_station_model.py
```

Supplementary figures are rebuilt with:

```bash
conda run -n japan-station-terms python build_event_adjusted_supplement_figures.py
```

The compact initial-submission files and peer-review code archive are built with:

```bash
conda run -n japan-station-terms python ../../work/build_cee_initial_submission_package.py
```

## Compile

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_information.tex
```

The primary MF2013 residual and hazard calculations use ground-surface installations and direct RotD100 coordinates. The Zhao 2006 sensitivity uses RotD50 as the closest available orientation-independent coordinate to its geometric-mean target. The ESM analysis uses post-2015 events within the published Bindi magnitude range and retains the failed Japanese-to-European transfer and adverse local spatial results. Hazard propagation uses zero-centred out-of-fold station predictions and converts the official `Vs=400 m/s` spectrum to the station AVS30 reference before applying the station term. All-earthquake, active-shallow and subduction products are evaluated separately; fixed-probability category ordinates are not summed. The resulting correction represents the source-path distribution sampled by the strong-motion archive. Empirical bounds describe station-model prediction error; frozen-network agreement is a separate reproducibility check. Neither is a full PSHA uncertainty interval.
