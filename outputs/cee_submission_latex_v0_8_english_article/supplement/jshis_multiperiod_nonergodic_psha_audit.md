# All-period official-response-spectrum non-ergodic PSHA audit

## Completed scope

- Station correction anchors: 9,044 station-target rows across PGA, SA(0.3 s), SA(1.0 s), and SA(3.0 s).
- Response-spectrum anchor stations: 2,261.
- Sampled-grid cells: 10,467.
- Official station response rows after correction: 72,352.
- Official sampled-grid response rows after correction: 334,944.
- Official response-spectrum periods: 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, and 5.0 s.

## Station validation

- SA(0.3s) RotD50: weighted RMSE 0.1533 log10, RMSE reduction 51.9% versus zero correction and 51.0% versus national mean.
- SA(1.0s) RotD50: weighted RMSE 0.0873 log10, RMSE reduction 61.9% versus zero correction and 50.1% versus national mean.
- SA(3.0s) RotD50: weighted RMSE 0.0678 log10, RMSE reduction 70.7% versus zero correction and 32.5% versus national mean.

## 50-year 10% sampled-grid multiplier spectrum

| period_s | source | q05 | q50 | q95 | delta_pct_q50 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0.1 | clamped_to_0.3s | 0.587 | 1.045 | 2.107 | 4.5 |
| 0.2 | clamped_to_0.3s | 0.587 | 1.045 | 2.107 | 4.5 |
| 0.3 | direct_anchor | 0.587 | 1.045 | 2.107 | 4.5 |
| 0.5 | log_period_interpolated | 0.593 | 0.966 | 1.553 | -3.4 |
| 1.0 | direct_anchor | 0.570 | 0.834 | 1.217 | -16.6 |
| 2.0 | log_period_interpolated | 0.545 | 0.720 | 0.992 | -28.0 |
| 3.0 | direct_anchor | 0.513 | 0.658 | 0.912 | -34.2 |
| 5.0 | clamped_to_3.0s | 0.513 | 0.658 | 0.912 | -34.2 |

## Boundary

This audit propagates non-ergodic station terms through official J-SHIS response-spectrum map ordinates at the matched 250 m meshes. The 0.3, 1.0, and 3.0 s corrections are direct residual-model anchors. The 0.1 and 0.2 s corrections are clamped to the 0.3 s anchor; the 0.5 and 2.0 s corrections are log-period interpolations; the 5.0 s correction is clamped to the 3.0 s anchor. This closes an all-period official-map sensitivity calculation. It does not close a source-level official J-SHIS PSHA rerun because the complete source recurrence, period-dependent sigma, spatial correlation, and logic-tree implementation are not present in the local public-data package.

## Output tables

- `outputs/jshis_multiperiod_station_corrections.csv`
- `outputs/jshis_multiperiod_station_surface_validation.csv`
- `outputs/jshis_multiperiod_station_surface_validation_summary.csv`
- `outputs/jshis_multiperiod_continuous_correction_surface_mesh.csv`
- `outputs/jshis_multiperiod_official_response_station_values.csv`
- `outputs/jshis_multiperiod_official_response_mesh_values.csv`
- `outputs/jshis_multiperiod_nonergodic_psha_summary.csv`
