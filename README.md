# Japan Strong Motion Station Terms

This repository contains the reproducible analysis and manuscript package for **Repeatable station terms modify fixed-probability long-period response spectra at stations across Japan**.

The study starts from 333,808 public K-NET and KiK-net records. The archive contains 231,380 ground-surface records; requiring finite F-net moment magnitude and positive AVS30 leaves 222,664 records at each of eight directly observed periods. The primary MF2013 calculation uses model-compatible RotD100 spectra. Residuals are decomposed into a global intercept, between-event residuals, site-to-site residuals (station terms) and within-site residuals. The station terms are tested across independent earthquake groups, a Zhao 2006 calculation, 102,428 paired KiK-net surface-to-borehole spectra, a fixed K-NET-to-KiK-net relation with disjoint earthquakes, spatial blocks and source/path strata. A post-2015 ESM v3 analysis supplies an independent European-Mediterranean sample. Source-conditioned terms separate crustal, interplate and intraplate records, and a publisher-supplied Kanto pattern provides an independent 3D-simulation comparison at 5.0 s. Cross-validated station terms are applied to matched J-SHIS response-spectrum products after conversion from the `Vs=400 m/s` engineering-bedrock reference to each station's AVS30 condition.

## Main result

At 3.0 s, station terms from separate two-way decompositions of independent Japanese earthquake groups correlate at `0.890`, and MF2013 and Zhao 2006 station terms correlate at `0.770`. A relation estimated from K-NET stations and four earthquake groups estimates KiK-net terms from the fifth group with correlation `0.598` and RMSE gain `20.8%`. The chronological K-NET-to-KiK-net result is `0.536` and `19.4%`. Paired surface-to-borehole response-spectrum-ratio and station-term shapes correlate at `0.635` when calculated from disjoint earthquakes. Source-matched Kanto 5.0 s station terms correlate at `0.304` with the independent simulation pattern; the 31.7 km block interval is `0.163--0.434`. Intraplate 3.0--5.0 s terms reduce spatial error by `14.3--14.8%`, while adverse source-class periods remain in the release. In 13,430 post-2015 ESM records, independent-earthquake-group station terms correlate at `0.928` at 3.0 s. Primary Japanese multipliers span `0.624--1.338` between the 5th and 95th percentiles.

## Repository layout

- `work/jshis_event_adjusted_station_model.py`: primary residual, validation and response-spectrum analysis.
- `work/jshis_mf2013_applicability_audit.py`: sensitivity to the complete published MF2013 regression-domain screens.
- `work/jshis_kiknet_surface_borehole_validation.py`: national paired-sensor spectral and independent-earthquake-group validation.
- `work/jshis_cross_network_transfer_validation.py`: cross-network and independent-earthquake station-term validation.
- `work/esm_external_station_validation.py`: post-2015 ESM recurrence, backbone sensitivity and adverse cross-region transfer.
- `work/jshis_source_category_response_sensitivity.py`: matched-cell propagation through official all-earthquake, active-shallow and subduction products.
- `work/jshis_source_conditioned_station_fields.py`: source-conditioned station terms and separate active-shallow/interplate/intraplate response-spectrum sensitivities.
- `work/sung2025_kanto_external_comparison.py`: source-matched SA(5 s) comparison with the independent Kanto 3D-simulation pattern of Sung et al. (2025).
- `work/jshis_hazard_impact_robustness.py`: interval calibration and independent station-term stability.
- `work/jshis_independent_gmpe_replication.py`: Zhao 2006 ground-motion-model sensitivity.
- `work/jshis_station_uncertainty_propagation.py`: empirical station-model intervals and spectrum propagation.
- `work/jshis_robustness_stress_tests.py`: network, macroregion and influential-event tests.
- `work/jshis_path_stratification_audit.py`: connected-component path and source stratification with event split-half checks.
- `work/jshis_balanced_station_model.py`: equal-stratum station terms and spatial-estimation sensitivity.
- `work/validate_event_adjusted_release.py`: release integrity checks used by CI.
- `environment.yml`: tested Python environment.
- `outputs/cee_submission_latex_v0_8_english_article/`: manuscript, Supplementary Information, figures and derived tables.

## Public inputs

Download the following official files and place them at these relative paths:

```text
work/external_data/jshis_gmf/flatfile_sub1-v2024.zip
work/external_data/jshis_mf2013/MF13rev_coefs.csv
work/external_data/jshis_respmap/P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip
work/external_data/jshis_respmap/P-Y2020-RESP-MAP-AVR-LND_MTTL-T50.zip
work/external_data/jshis_respmap/P-Y2020-RESP-MAP-AVR-PPE_MTTL-T50.zip
work/external_data/sung2025/bssa-2024239_supplement.csv.zip
```

Sources:

- J-SHIS/NIED Strong Ground Motion Flat File: <https://doi.org/10.17598/NIED.0032>
- MF2013 coefficients: <https://www.j-shis.bosai.go.jp/labs/mf2013/en/>
- J-SHIS response-spectrum products: <https://www.j-shis.bosai.go.jp/map/respmap/>
- Sung et al. (2025) Kanto T=5 s adjusted site field: <https://doi.org/10.1785/0120240239>

Raw third-party data are not redistributed.

Extract `BSSA-2024239_Supplement.csv` from the Sung et al. archive into the same `work/external_data/sung2025/` directory before running the external-field comparison.

The ESM workflow downloads its exact version 3 flatfile query to the ignored path `work/external_data/esm_v3/esm_v3_external_validation.csv`. The analysed snapshot has SHA-256 `084ef4c7519fdeb26a05db9ad980535cd4c2207b2358d6b668cd5dafb4f0c724`; the script stores the query URL and digest with the local input metadata. ESM v3 is available at <https://doi.org/10.13127/esm.3>.

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
conda run -n japan-station-terms python work/esm_external_station_validation.py
conda run -n japan-station-terms python work/jshis_source_category_response_sensitivity.py
conda run -n japan-station-terms python work/jshis_source_conditioned_station_fields.py
conda run -n japan-station-terms python work/sung2025_kanto_external_comparison.py
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

The response-spectrum calculation is a matched-surface-station sensitivity analysis for the source-path distribution sampled by the strong-motion archive. It retains the source recurrence, aleatory model and logic-tree aggregation within each published J-SHIS product. Fixed-probability active-shallow and subduction ordinates are evaluated separately and are not summed. Source-conditioned crustal, interplate and intraplate fields are separate sensitivities; they are not combined without official source weights. The MF2013 applicability sensitivity reconstructs the published Kanno-PGA screening equations on the current public flatfile and does not recreate the historical waveform database or regression weights. The adverse ESM transfer, local spatial and source-conditioned results are retained. Empirical bounds quantify station-model prediction error, while frozen-network agreement measures reproducibility under a disjoint network and earthquake set. The release does not claim an official source-level J-SHIS PSHA reproduction.
