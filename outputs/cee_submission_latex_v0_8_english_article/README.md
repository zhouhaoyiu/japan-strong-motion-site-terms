# CEE Manuscript Package

## Manuscript

- `main.tex` and `main.pdf`: CEE Article manuscript.
- `supplementary_information.tex` and `supplementary_information.pdf`: Supplementary Information.
- `cover_letter_cee.md`: cover-letter draft.
- `advisor_summary_cn.md`: Chinese project summary for advisor review.

## Current figures

The main manuscript uses six figures:

1. `figure_event_adjusted_overview`
2. `figure_event_adjusted_site_structure`
3. `figure_ground_motion_model_sensitivity`
4. `figure_event_adjusted_multiperiod`
5. `figure_path_stratification`
6. `figure_event_adjusted_city_cases`

The Supplementary Information uses seven figures: event repeatability, spatial validation, robustness stress tests, equal-stratum prediction, response-spectrum sensitivity, station-model uncertainty and the KiK-net waveform comparison. Vector PDFs are the manuscript sources; PNG copies are retained for portal compatibility.

## Derived tables

The `supplement/` directory contains the primary MF2013 results, the RotD50--RotD100 component audit, Zhao 2006 model sensitivity, station-model prediction intervals, transfer tests, influential-event deletion results, path/source stratification and equal-stratum sensitivity. Rebuild the analyses in this order:

```bash
conda run -n japan-station-terms python ../../work/jshis_event_adjusted_station_model.py
conda run -n japan-station-terms python ../../work/jshis_independent_gmpe_replication.py
conda run -n japan-station-terms python ../../work/jshis_station_uncertainty_propagation.py
conda run -n japan-station-terms python ../../work/jshis_robustness_stress_tests.py
conda run -n japan-station-terms python ../../work/jshis_path_stratification_audit.py
conda run -n japan-station-terms python ../../work/jshis_balanced_station_model.py
```

Supplementary figures are rebuilt with:

```bash
conda run -n japan-station-terms python build_event_adjusted_supplement_figures.py
```

## Compile

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_information.tex
```

The primary MF2013 residual and hazard calculations use ground-surface installations and direct RotD100 coordinates. The Zhao 2006 sensitivity uses RotD50 as the closest available orientation-independent coordinate to its geometric-mean target. Hazard propagation uses zero-centred out-of-fold station predictions and converts the official `Vs=400 m/s` spectrum to the station AVS30 reference before applying the station term. The resulting correction represents the source-path distribution sampled by the strong-motion archive. Empirical bounds describe station-model prediction error and are not full PSHA uncertainty intervals.
