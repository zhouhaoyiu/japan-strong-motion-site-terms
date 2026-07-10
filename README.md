# Japan Strong Motion Station Terms

This repository contains the reproducible analysis and manuscript package for **Repeatable station terms redistribute long-period response spectra across Japanese strong-motion sites**.

The study starts from 333,808 public K-NET and KiK-net records and retains 222,664 ground-surface records at each of eight directly observed periods. The primary MF2013 calculation uses model-compatible RotD100 spectra. Residuals are decomposed into global, event and station components, then tested across held-out events, a Zhao 2006 calculation, spatial blocks, networks, fixed macroregions, influential-event deletions and path/source strata. Cross-validated station predictions and empirical prediction intervals are applied to matched J-SHIS response-spectrum ordinates after converting the official `Vs=400 m/s` reference to each station's AVS30 condition.

## Main result

At 3.0 s, station terms estimated from disjoint event groups correlate at `0.921`, and MF2013 and Zhao 2006 station terms correlate at `0.770`. Site and location variables reduce spatial-block RMSE by `12.6%` across 1,628 surface stations. Directional-sector correlations with the full station field range from `0.280` to `0.886`, showing path-conditioned transfer. Cross-validated multipliers have 5th and 95th percentiles of `0.624` and `1.338`. The matched-station median 50-year 10% surface spectrum changes from `0.081 g` to `0.075 g`; the nominal 90% empirical station-model interval attains `89.0%` station coverage.

## Repository layout

- `work/jshis_event_adjusted_station_model.py`: primary residual, validation and response-spectrum analysis.
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

## Reproduction

```bash
conda env create -f environment.yml
conda run -n japan-station-terms python work/jshis_event_adjusted_station_model.py
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

The response-spectrum calculation is a matched-surface-station sensitivity analysis for the source-path distribution sampled by the strong-motion archive. It retains the source recurrence, aleatory model and logic-tree aggregation of the published J-SHIS product. The empirical bounds quantify station-model prediction error only. The release does not claim an official source-level J-SHIS PSHA reproduction, source-specific transfer or a complete PSHA uncertainty interval.
