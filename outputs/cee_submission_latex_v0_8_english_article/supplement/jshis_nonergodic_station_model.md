# Non-ergodic Station Residual Model Experiment

Status: station-held-out supplement for the J-SHIS/MF2013 residual manuscript.

## Design

- Input station residuals: `mf2013_site` station mean residuals with at least 20 records per target and station.
- Target: station mean residual after the official MF2013 D1400/AVS30 site correction.
- Baseline: zero additional station correction.
- Models: ridge regression with public site variables, ridge regression with site plus coordinates, and gradient-boosted regression with site plus coordinates.
- Validation: five-fold random station holdout and five-fold spatial-block holdout based on station longitude and latitude.
- Metrics: weighted station-residual MAE/RMSE, weighted by the number of records per station.

## Best model by target

### Random station holdout

| target_label | model | rmse_reduction_pct_mean | rmse_reduction_pct_std | mae_reduction_pct_mean | weighted_r2_vs_zero_mean |
| --- | --- | --- | --- | --- | --- |
| PGA RotD50 | gradient_boosted_site_space | 60.912 | 2.010 | 63.602 | 0.847 |
| SA(0.3s) RotD50 | gradient_boosted_site_space | 54.124 | 3.027 | 56.013 | 0.789 |
| SA(1.0s) RotD50 | gradient_boosted_site_space | 53.496 | 4.154 | 56.013 | 0.782 |
| SA(3.0s) RotD50 | gradient_boosted_site_space | 61.350 | 3.483 | 65.610 | 0.850 |

### Spatial-block holdout

| target_label | model | rmse_reduction_pct_mean | rmse_reduction_pct_std | mae_reduction_pct_mean | weighted_r2_vs_zero_mean |
| --- | --- | --- | --- | --- | --- |
| PGA RotD50 | gradient_boosted_site_space | 45.622 | 7.496 | 47.266 | 0.700 |
| SA(0.3s) RotD50 | gradient_boosted_site_space | 41.741 | 4.329 | 42.363 | 0.659 |
| SA(1.0s) RotD50 | gradient_boosted_site_space | 38.898 | 8.007 | 41.585 | 0.622 |
| SA(3.0s) RotD50 | gradient_boosted_site_space | 44.345 | 8.960 | 48.633 | 0.684 |

## Interpretation

The spatial-block test retains positive long-period residual reduction, so the supplement supports a predictive non-ergodic site-residual extension beyond the official D1400/AVS30 correction.

This experiment does not define a new ground-motion prediction equation. It tests whether public site and location variables can predict remaining station terms under held-out-station validation.

## Output files

- `jshis_nonergodic_station_model_fold_metrics.csv`: fold-level validation metrics.
- `jshis_nonergodic_station_model_summary.csv`: summary metrics by target, split, and model.
- `jshis_nonergodic_station_model_predictions.csv`: held-out station predictions.
- `figures/cee_fig6_nonergodic_station_model.pdf`: manuscript-style summary figure.
