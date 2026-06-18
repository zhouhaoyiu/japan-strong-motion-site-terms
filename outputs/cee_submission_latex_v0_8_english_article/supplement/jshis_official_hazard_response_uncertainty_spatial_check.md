# Official J-SHIS Response-Spectrum Uncertainty and Spatial Check

- Matched stations at 50y_10pct: 2261
- Station multiplier q05/q50/q95: 0.414 / 0.619 / 0.862
- Official SA(3.0 s) q05/q50/q95: 0.034 / 0.076 / 0.202 g
- Corrected SA(3.0 s) q05/q50/q95: 0.019 / 0.048 / 0.130 g
- Spatial-block corrected median range across folds: 0.037--0.061 g

This check uses station-level multipliers from the spatial-block held-out residual model. It adjusts only the SA(3.0 s) ordinate from the official J-SHIS response-spectrum map. It is not a full PSHA rerun.

## Bootstrap confidence intervals

| probability_level | metric | distribution_quantile | observed | bootstrap_ci_low | bootstrap_ci_high | n_stations | bootstrap_reps | seed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50y_10pct | station_multiplier | 0.05 | 0.413627 | 0.40371 | 0.426818 | 2261 | 2000 | 20260612 |
| 50y_10pct | station_multiplier | 0.5 | 0.619055 | 0.61472 | 0.625408 | 2261 | 2000 | 20260612 |
| 50y_10pct | station_multiplier | 0.95 | 0.861878 | 0.845882 | 0.875234 | 2261 | 2000 | 20260612 |
| 50y_10pct | official_sa3_g | 0.05 | 0.0338561 | 0.0333846 | 0.0349549 | 2261 | 2000 | 20260612 |
| 50y_10pct | official_sa3_g | 0.5 | 0.0755655 | 0.0733382 | 0.0779444 | 2261 | 2000 | 20260612 |
| 50y_10pct | official_sa3_g | 0.95 | 0.20228 | 0.192001 | 0.210264 | 2261 | 2000 | 20260612 |
| 50y_10pct | corrected_sa3_g | 0.05 | 0.0192505 | 0.0186973 | 0.0200596 | 2261 | 2000 | 20260612 |
| 50y_10pct | corrected_sa3_g | 0.5 | 0.0475316 | 0.0463023 | 0.0493539 | 2261 | 2000 | 20260612 |
| 50y_10pct | corrected_sa3_g | 0.95 | 0.12977 | 0.124691 | 0.135543 | 2261 | 2000 | 20260612 |
| 50y_10pct | delta_pct | 0.05 | -58.6373 | -59.629 | -57.4144 | 2261 | 2000 | 20260612 |
| 50y_10pct | delta_pct | 0.5 | -38.0945 | -38.5462 | -37.4194 | 2261 | 2000 | 20260612 |
| 50y_10pct | delta_pct | 0.95 | -13.8122 | -15.55 | -12.4382 | 2261 | 2000 | 20260612 |
