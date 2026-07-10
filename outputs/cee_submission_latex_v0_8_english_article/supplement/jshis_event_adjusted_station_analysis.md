# Event-adjusted direct-period station analysis

## Primary design

- RotD100 response spectra at eight periods (0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, and 5.0 s) are read directly from the public flatfile to match the MF2013 horizontal-vector response definition.
- Each residual is decomposed into a global intercept, a zero-centred event term, a zero-centred station term, and a record remainder.
- Station prediction uses five spatial blocks and only out-of-fold predictions in the hazard sensitivity calculation.
- Official response-map ordinates are retained at Vs=400 m/s and are separately projected to each station AVS30 before the station adjustment is applied.

## SA(3.0 s) result

- Supported stations: 1,628.
- Spatial-block RMSE reduction versus a zero station term: 12.6%.
- Observed-predicted station-term correlation: 0.500.
- Across held-out event groups, the mean train-test station-term correlation is 0.921 and the mean test RMSE reduction is 74.5%.
- Cross-validated station multiplier q05/q50/q95: 0.624 / 0.886 / 1.338.
- Matched-station 50-year 10% official Vs400 median: 0.078 g.
- MF2013 surface-reference median: 0.081 g.
- Station-adjusted surface-reference median: 0.075 g.

## Model boundary

- The primary attenuation backbone contains the MF2013 basic, D1400, and AVS30 terms. It is labelled as such and is not described as the complete official implementation.
- The rule-based AI sensitivity gives an SA(3.0 s) station-term correlation of 0.989 with the primary decomposition.
- PH is an event-constant period term for qualifying Philippine Sea Plate intraplate earthquakes. The event fixed effect absorbs it for station-term estimation; the global intercept is not interpreted as an official MF2013 bias.
- The response-map calculation is a matched-station surface-spectrum sensitivity analysis, not an official source-level J-SHIS PSHA rerun.

## Output files

- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_adjusted_station_terms.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_station_decomposition_metrics.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_adjusted_station_model_predictions.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_adjusted_station_model_metrics.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_holdout_station_repeatability.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_adjusted_surface_spectrum_values.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_adjusted_surface_spectrum_summary.csv`
- `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_event_adjusted_city_nearest_cases.csv`
