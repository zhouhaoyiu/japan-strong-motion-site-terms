# Additional CEE submission experiments

Generated on 2026-06-14. These tables cover the required, recommended, and feasible optional checks requested before CEE submission.

## Regional extrapolation

| target_label | period_label | model | split | folds | n_test_sites_total | n_test_sites_mean | weighted_mae_mean | weighted_mae_sd | weighted_rmse_mean | weighted_rmse_sd | mae_reduction_pct_mean | mae_reduction_pct_min | mae_reduction_pct_max | rmse_reduction_pct_mean | rmse_reduction_pct_min | rmse_reduction_pct_max | weighted_r2_vs_zero_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SA(3.0s) RotD50 | SA(3.0 s) | full_site_location | region_leaveout | 7 | 2261 | 323 | 0.1145 | 0.01512 | 0.1478 | 0.01918 | 52.36 | 41.55 | 64.94 | 47.82 | 37.9 | 59.88 | 0.722 |

## Official-variable substitution for SA(3.0 s)

| target_label | period_label | model | split | folds | n_test_sites_total | n_test_sites_mean | weighted_mae_mean | weighted_mae_sd | weighted_rmse_mean | weighted_rmse_sd | mae_reduction_pct_mean | mae_reduction_pct_min | mae_reduction_pct_max | rmse_reduction_pct_mean | rmse_reduction_pct_min | rmse_reduction_pct_max | weighted_r2_vs_zero_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SA(3.0s) RotD50 | SA(3.0 s) | full_site_location | spatial_block | 5 | 2261 | 452.2 | 0.116 | 0.01668 | 0.1506 | 0.02043 | 51 | 41.32 | 60.57 | 46.25 | 36.1 | 55.36 | 0.7056 |
| SA(3.0s) RotD50 | SA(3.0 s) | mf2013_d1400_avs30 | spatial_block | 5 | 2261 | 452.2 | 0.1423 | 0.01855 | 0.189 | 0.02241 | 39.97 | 30 | 48.48 | 32.75 | 23.36 | 39.79 | 0.5438 |
| SA(3.0s) RotD50 | SA(3.0 s) | official_four | spatial_block | 5 | 2261 | 452.2 | 0.1308 | 0.01701 | 0.1738 | 0.02551 | 44.82 | 33.36 | 53.28 | 38.24 | 27.47 | 46.61 | 0.6138 |

## Positive SA(3.0 s) corrections

| metric | value | unit | note |
| --- | --- | --- | --- |
| all_stations | 2261 | stations | Stations with held-out SA(3.0 s) correction and official response-map match. |
| positive_correction_stations | 24 | stations | Station multiplier greater than 1.0 at 50-year 10% probability. |
| positive_correction_share_pct | 1.061 | percent | Share of matched stations with upward correction. |
| positive_multiplier_median | 1.064 | ratio | Median upward multiplier among positive-correction stations. |
| positive_multiplier_95pct | 1.333 | ratio | 95th percentile upward multiplier among positive-correction stations. |
| positive_multiplier_max | 1.55 | ratio | Largest upward multiplier among positive-correction stations. |
| positive_share_Chubu | 1.899 | percent | 9/474 matched stations in Chubu. |
| positive_share_Chugoku_Shikoku | 0.3378 | percent | 1/296 matched stations in Chugoku_Shikoku. |
| positive_share_Hokkaido | 0.303 | percent | 1/330 matched stations in Hokkaido. |
| positive_share_Kanto | 0.3425 | percent | 1/292 matched stations in Kanto. |
| positive_share_Kinki | 0 | percent | 0/189 matched stations in Kinki. |
| positive_share_Kyushu_Okinawa | 0.3509 | percent | 1/285 matched stations in Kyushu_Okinawa. |

## Matching audit

| check | count | denominator | percent | source | note |
| --- | --- | --- | --- | --- | --- |
| strong_motion_records | 333808 |  |  | jshis_smrec_schema_sub1_summary.md | Public flatfile records in sub1-v2024. |
| unique_record_sites | 2581 |  |  | jshis_smrec_schema_sub1_summary.md | Unique site_id//10 values in the strong-motion records. |
| site_schema_rows | 2607 |  |  | jshis_site_schema_v2024_sub1.csv | Rows in the public site schema. |
| record_site_schema_matched_rows | 333808 | 3.338e+05 | 100 | jshis_smrec_schema_sub1_summary.md | Record-to-site match reported by the schema extraction audit. |

Full CSV outputs contain all periods and all feature sets.
