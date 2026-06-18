# SA(3.0 s) physical-proxy station model

This check fits a transparent ridge model to the MF2013 site-corrected SA(3.0 s) station residuals. Predictors are basin-depth proxies, shallow-velocity proxies, public borehole-profile proxies, volcanic-front distances, elevation, sensor depth, coordinates, and broad geographic region. No station identifier enters the model.

## Spatial-block held-out performance

| target_label | n_stations | weighted_rmse_zero | weighted_rmse_national_mean | weighted_rmse_physical_proxy | rmse_reduction_vs_zero_pct | rmse_reduction_vs_mean_pct | pearson_observed_predicted |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SA(3.0s) RotD50 | 2261 | 0.2931 | 0.1922 | 0.1915 | 34.6757 | 0.3763 | 0.2935 |

## Largest standardized coefficients

| feature | standardized_coefficient |
| --- | --- |
| region_Tohoku | 0.0659 |
| log_vs20 | -0.0643 |
| log_dbase | 0.0576 |
| region_Kyushu_Okinawa | -0.056 |
| lat | 0.0412 |
| lon | -0.0379 |
| log_avs30 | 0.037 |
| t_d1400_vs700_s | -0.0353 |
| region_Hokkaido | -0.0312 |
| depth_to_vs_ge_500_m | -0.0152 |
| dist_vf_mf13_nejapan | -0.0125 |
| region_Kanto | 0.0106 |

## Interpretation boundary

The model is an independent physical and geographic proxy check. Calibrated three-dimensional basin-response simulation requires separate velocity and basin-geometry data. The proxy model shows that a simple interpretable model recovers part of the long-period station term outside the gradient-boosted station model.
