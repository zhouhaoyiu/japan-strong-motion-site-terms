# Official J-SHIS Hazard and Response-Spectrum Check

This check applies the held-out SA(3.0 s) station correction to official J-SHIS 2020 response-spectrum map ordinates at the same 250 m station mesh. The response-spectrum data are on engineering bedrock (Vs=400 m/s) and use 5% damping. The calculation adjusts only the SA(3.0 s) ordinate and does not constitute a full PSHA rerun.

## Data sources

- J-SHIS PSHM hazard-curve API, Y2024, average case, all earthquakes, 30-year and 50-year periods.
- J-SHIS response-spectrum map data, Y2020, average case, all earthquakes, 50-year period.

## Station-level official SA(3.0 s) map values

- Matched station-probability rows: 9044
- Matched held-out stations at 50-year 10% level: 2261
- 50-year 10% official SA(3.0 s) q05/q50/q95: 0.034 / 0.076 / 0.202 g
- 50-year 10% corrected SA(3.0 s) q05/q50/q95: 0.019 / 0.048 / 0.130 g

## Example stations at 50-year 10%

| example_type | site_code | meshcode250 | official_SA3_g | corrected_SA3_g | multiplier |
| --- | --- | --- | --- | --- | --- |
| maximum_negative | KYTH07 | 5235257943 | 0.144 | 0.023 | 0.163 |
| q05 | KOCH10 | 5033233533 | 0.150 | 0.062 | 0.414 |
| q25 | MIE012 | 5136360542 | 0.290 | 0.155 | 0.537 |
| q50 | FKS031 | 5640060423 | 0.073 | 0.045 | 0.619 |
| q75 | FKO008 | 5030268832 | 0.049 | 0.035 | 0.706 |
| q95 | MYZH15 | 4831443733 | 0.075 | 0.065 | 0.862 |
| maximum_positive | NIG024 | 5538535523 | 0.051 | 0.079 | 1.550 |

## Summary table

| metric | n_station_values | official_sa3_g_q05 | official_sa3_g_q50 | official_sa3_g_q95 | corrected_sa3_g_q05 | corrected_sa3_g_q50 | corrected_sa3_g_q95 | delta_pct_q05 | delta_pct_q50 | delta_pct_q95 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| official_sa3_50y_2pct | 2261 | 0.053 | 0.119 | 0.304 | 0.030 | 0.074 | 0.202 | -58.637 | -38.094 | -13.812 |
| official_sa3_50y_5pct | 2261 | 0.042 | 0.092 | 0.241 | 0.024 | 0.058 | 0.156 | -58.637 | -38.094 | -13.812 |
| official_sa3_50y_10pct | 2261 | 0.034 | 0.076 | 0.202 | 0.019 | 0.048 | 0.130 | -58.637 | -38.094 | -13.812 |
| official_sa3_50y_39pct | 2261 | 0.021 | 0.045 | 0.125 | 0.012 | 0.028 | 0.080 | -58.637 | -38.094 | -13.812 |
| example_station_multipliers | 7 |  |  |  |  |  |  | -76.213 | -38.094 | 34.341 |
| official_api_curve_points | 264 |  |  |  |  |  |  |  |  |  |
