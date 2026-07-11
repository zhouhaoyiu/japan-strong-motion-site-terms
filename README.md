# Japan Strong Motion Station Terms

This repository contains the reproducible analysis and manuscript package for **Repeatable station terms redistribute long-period response spectra across Japanese strong-motion sites**.

The study starts from 333,808 public K-NET and KiK-net records. The archive contains 231,380 ground-surface records; requiring finite F-net moment magnitude and positive AVS30 leaves 222,664 records at each of eight directly observed periods. The primary MF2013 calculation uses model-compatible RotD100 spectra. Residuals are decomposed into global, event and station components, then tested across held-out events, a Zhao 2006 calculation, 102,428 paired KiK-net surface--borehole spectra, a frozen K-NET-to-KiK-net prediction with disjoint earthquakes, spatial blocks, fixed macroregions, influential-event deletions and path/source strata. Cross-validated station predictions and empirical prediction intervals are applied to matched J-SHIS response-spectrum ordinates after converting the official `Vs=400 m/s` reference to each station's AVS30 condition.

## Main result

At 3.0 s, station terms from separate two-way decompositions of disjoint event groups correlate at `0.890`, and MF2013 and Zhao 2006 station terms correlate at `0.770`. A model trained on K-NET stations and training earthquakes predicts KiK-net terms from held-out earthquakes with correlation `0.598` and RMSE gain `20.8%`. When events are ordered by origin time, the earliest 80% of K-NET events predict KiK-net terms from the latest 20% with correlation `0.536` and RMSE gain `19.4%`; all eight period-specific lower confidence bounds remain positive. Paired surface--borehole and station-term spectral shapes correlate at `0.635` when estimated from disjoint events. The frozen K-NET and primary correction fields correlate at `0.736` and agree in direction at `79.7%` of KiK-net sites. Primary cross-validated multipliers have 5th and 95th percentiles of `0.624` and `1.338`; a calibrated 90% prediction interval excludes unity at `2.0%` of stations. The complete MF2013 regression-domain screen retains 35,857 records and a station-field correlation of `0.959`.

## Repository layout

- `work/jshis_event_adjusted_station_model.py`: primary residual, validation and response-spectrum analysis.
- `work/jshis_mf2013_applicability_audit.py`: sensitivity to the complete published MF2013 regression-domain screens.
- `work/jshis_kiknet_surface_borehole_validation.py`: national paired-sensor spectral and held-event validation.
- `work/jshis_cross_network_transfer_validation.py`: network- and event-external frozen-model validation.
- `work/jshis_hazard_impact_robustness.py`: interval calibration and independent correction-field stability.
- `work/jshis_independent_gmpe_replication.py`: Zhao 2006 ground-motion-model sensitivity.
- `work/jshis_station_uncertainty_propagation.py`: empirical station-model intervals and spectrum propagation.
- `work/jshis_robustness_stress_tests.py`: network, macroregion and influential-event tests.
- `work/jshis_path_stratification_audit.py`: connected-component path and source stratification with event split-half checks.
- `work/jshis_balanced_station_model.py`: equal-stratum station fields and spatial prediction sensitivity.
- `work/validate_event_adjusted_release.py`: release integrity checks used by CI.
- `environment.yml`: tested Python environment.
- `outputs/cee_submission_latex_v0_8_english_article/`: manuscript, Supplementary Information, figures and derived tables.

## Public inputs

Download the following official files and place them at these relative paths:

```text
work/external_data/jshis_gmf/flatfile_sub1-v2024.zip
work/external_data/jshis_mf2013/MF13rev_coefs.csv
work/external_data/jshis_respmap/P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip
```

Sources:

- J-SHIS/NIED Strong Ground Motion Flat File: <https://doi.org/10.17598/NIED.0032>
- MF2013 coefficients: <https://www.j-shis.bosai.go.jp/labs/mf2013/en/>
- J-SHIS response-spectrum products: <https://www.j-shis.bosai.go.jp/>

Raw third-party data are not redistributed.

`public_inputs_manifest.tsv` records the exact file sizes, SHA-256 hashes, source pages and required archive members used for this release. Verify downloaded inputs before running the analysis:

```bash
conda run -n japan-station-terms python work/verify_public_inputs.py
```

## Reproduction

```bash
conda env create -f environment.yml
conda run -n japan-station-terms python work/verify_public_inputs.py
conda run -n japan-station-terms python work/audit_jshis_flatfile_selection.py
conda run -n japan-station-terms python work/jshis_event_adjusted_station_model.py
conda run -n japan-station-terms python work/jshis_spatial_model_complexity_audit.py
conda run -n japan-station-terms python work/jshis_mf2013_applicability_audit.py
conda run -n japan-station-terms python work/jshis_kiknet_surface_borehole_validation.py
conda run -n japan-station-terms python work/jshis_cross_network_transfer_validation.py
conda run -n japan-station-terms python work/jshis_hazard_impact_robustness.py
conda run -n japan-station-terms python work/jshis_independent_gmpe_replication.py
conda run -n japan-station-terms python work/jshis_station_uncertainty_propagation.py
conda run -n japan-station-terms python work/jshis_robustness_stress_tests.py
conda run -n japan-station-terms python work/jshis_path_stratification_audit.py
conda run -n japan-station-terms python work/jshis_balanced_station_model.py
conda run -n japan-station-terms python outputs/cee_submission_latex_v0_8_english_article/build_event_adjusted_supplement_figures.py
python work/validate_event_adjusted_release.py
```

Compile the manuscripts from the article directory:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_information.tex
```

## Scope

The response-spectrum calculation is a matched-surface-station sensitivity analysis for the source-path distribution sampled by the strong-motion archive. It retains the source recurrence, aleatory model and logic-tree aggregation of the published J-SHIS product. The MF2013 applicability sensitivity reconstructs the published Kanno-PGA screening equations on the current public flatfile; it does not recreate the historical waveform database or regression weights. Empirical bounds quantify station-model prediction error, while frozen-network agreement measures reproducibility under a disjoint network and earthquake set. Neither is a complete PSHA uncertainty interval. The release does not claim an official source-level J-SHIS PSHA reproduction or source-specific transfer.
