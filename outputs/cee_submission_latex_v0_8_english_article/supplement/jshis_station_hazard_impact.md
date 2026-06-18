# SA(3.0 s) Hazard-Curve and UHS Impact

This analysis is a hazard-impact transformation, not a full PSHA run. The local inputs do not include source recurrence rates or official site-specific hazard curves.

If the station correction is c in log10 spectral acceleration, the multiplicative SA(3.0 s) factor is f = 10^c. For any baseline hazard curve lambda0(S), the corrected curve is lambda_c(S) = lambda0(S / f). At a fixed annual exceedance rate, the SA(3.0 s) UHS ordinate becomes f times the baseline ordinate. For a local power-law hazard slope lambda0(S) proportional to S^-beta, the fixed-ordinate annual exceedance-rate ratio is f^beta.

## Main station-level distribution

- Stations: 2261
- SA(3.0 s) UHS multiplier q05/q50/q95: 0.414 / 0.619 / 0.862
- SA(3.0 s) UHS multiplier range: 0.163 to 1.550
- Fixed-SA hazard-rate ratio with beta=3 q05/q50/q95: 0.071 / 0.237 / 0.640

## Example stations

| example_type | site_code | network_label | correction_log10 | UHS_multiplier | hazard_ratio_beta3 |
| --- | --- | --- | --- | --- | --- |
| q05 | KOCH10 | KiK-net | -0.383 | 0.414 | 0.071 |
| q25 | MIE012 | K-NET | -0.270 | 0.537 | 0.155 |
| q50 | FKS031 | K-NET | -0.208 | 0.619 | 0.237 |
| q75 | FKO008 | K-NET | -0.151 | 0.706 | 0.352 |
| q95 | MYZH15 | KiK-net | -0.065 | 0.862 | 0.640 |
| maximum_positive | NIG024 | K-NET | 0.190 | 1.550 | 3.722 |
| maximum_negative | KYTH07 | KiK-net | -0.789 | 0.163 | 0.004 |
