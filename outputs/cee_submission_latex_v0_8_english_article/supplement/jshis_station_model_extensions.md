# Station Residual Model Extensions

## Ablation summary

| target_label | feature_set | rmse_reduction_pct_mean | rmse_reduction_pct_std | mae_reduction_pct_mean |
| --- | --- | --- | --- | --- |
| SA(1.0s) RotD50 | d1400_avs30 | 13.08 | 9.42 | 16.76 |
| SA(1.0s) RotD50 | depth_only | 11.80 | 10.53 | 15.07 |
| SA(1.0s) RotD50 | full_site_location | 38.26 | 7.97 | 40.35 |
| SA(1.0s) RotD50 | location_only | 13.21 | 10.94 | 16.15 |
| SA(1.0s) RotD50 | vs_only | 21.09 | 10.13 | 24.04 |
| SA(3.0s) RotD50 | d1400_avs30 | 33.42 | 6.96 | 40.40 |
| SA(3.0s) RotD50 | depth_only | 33.81 | 5.31 | 40.87 |
| SA(3.0s) RotD50 | full_site_location | 47.06 | 7.60 | 51.75 |
| SA(3.0s) RotD50 | location_only | 32.92 | 9.83 | 39.14 |
| SA(3.0s) RotD50 | vs_only | 38.76 | 7.82 | 45.57 |

## Grouped permutation importance

| target_label | group | rmse_increase_pct_mean | rmse_increase_pct_std |
| --- | --- | --- | --- |
| PGA RotD50 | station_metadata | 104.60 | 3.72 |
| PGA RotD50 | shallow_velocity | 13.94 | 2.21 |
| PGA RotD50 | volcanic_front_distance | 7.26 | 8.14 |
| PGA RotD50 | basin_depth | 2.45 | 2.75 |
| PGA RotD50 | location | -0.05 | 0.52 |
| SA(0.3s) RotD50 | station_metadata | 91.26 | 8.21 |
| SA(0.3s) RotD50 | shallow_velocity | 26.50 | 5.80 |
| SA(0.3s) RotD50 | volcanic_front_distance | 4.30 | 5.38 |
| SA(0.3s) RotD50 | basin_depth | 1.92 | 1.73 |
| SA(0.3s) RotD50 | location | 0.11 | 0.79 |
| SA(1.0s) RotD50 | station_metadata | 53.97 | 14.68 |
| SA(1.0s) RotD50 | shallow_velocity | 25.89 | 5.69 |
| SA(1.0s) RotD50 | basin_depth | 3.38 | 3.01 |
| SA(1.0s) RotD50 | volcanic_front_distance | 0.99 | 0.84 |
| SA(1.0s) RotD50 | location | -0.09 | 0.10 |
| SA(3.0s) RotD50 | station_metadata | 22.40 | 12.63 |
| SA(3.0s) RotD50 | shallow_velocity | 13.86 | 4.44 |
| SA(3.0s) RotD50 | basin_depth | 5.97 | 2.76 |
| SA(3.0s) RotD50 | volcanic_front_distance | 0.81 | 1.48 |
| SA(3.0s) RotD50 | location | -0.37 | 1.84 |

## K-NET/KiK-net transfer

| target_label | train_network | test_network | n_train_sites | n_test_sites | rmse_reduction_pct |
| --- | --- | --- | --- | --- | --- |
| SA(1.0s) RotD50 | K-NET | KiK-net | 987 | 1274 | 11.50 |
| SA(1.0s) RotD50 | KiK-net | K-NET | 1274 | 987 | 11.34 |
| SA(3.0s) RotD50 | K-NET | KiK-net | 987 | 1274 | 36.72 |
| SA(3.0s) RotD50 | KiK-net | K-NET | 1274 | 987 | 19.14 |

## Engineering-scale examples for SA(3.0 s)

| site_code | network_label | n_records | heldout_station_residual_log10 | predicted_station_correction_log10 | predicted_factor | residual_after_correction_log10 |
| --- | --- | --- | --- | --- | --- | --- |
| NIG024 | K-NET | 68 | 0.184 | 0.190 | 1.550 | -0.006 |
| AKTH11 | KiK-net | 26 | 0.269 | 0.130 | 1.350 | 0.139 |
| AKTH14 | KiK-net | 335 | -0.083 | 0.092 | 1.235 | -0.175 |
| YMN008 | K-NET | 45 | -0.288 | 0.081 | 1.205 | -0.369 |
| NIGH08 | KiK-net | 383 | 0.427 | 0.080 | 1.202 | 0.347 |
| AOM013 | K-NET | 222 | -0.130 | 0.058 | 1.142 | -0.187 |
| KYTH07 | KiK-net | 33 | -0.858 | -0.789 | 0.163 | -0.069 |
| OSKH01 | KiK-net | 32 | -1.006 | -0.773 | 0.169 | -0.233 |
| MYGH01 | KiK-net | 143 | -0.689 | -0.747 | 0.179 | 0.058 |
| OSKH02 | KiK-net | 31 | -1.086 | -0.732 | 0.185 | -0.353 |
| KYTH08 | KiK-net | 21 | -0.687 | -0.674 | 0.212 | -0.013 |
| IBUH06 | KiK-net | 125 | -0.731 | -0.652 | 0.223 | -0.079 |

## Interpretation

The ablation separates basin-depth, shallow-velocity, location, and full site-plus-location information. The grouped permutation analysis evaluates which predictor groups carry the held-out predictive signal. The network-transfer test is stricter than random station splits because training and testing use different Japanese strong-motion networks. The engineering examples translate log10 station corrections into multiplicative spectral-amplitude factors.
