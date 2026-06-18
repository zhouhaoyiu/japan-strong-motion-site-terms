# Sigma sensitivity audit for the 1,797-mesh national stratified recomputation

This audit compares five log10 sigma values on the same deterministic national stratified sample. The target is model calibration against official J-SHIS map probabilities at the official PBV thresholds, using the public-parameter recomputation pipeline.

## Integrity

- sigma=0.18: meshes=1,797, rows=7,188/7,188, failed=0, zero_probability=0, duplicates=0
- sigma=0.20: meshes=1,797, rows=7,188/7,188, failed=0, zero_probability=0, duplicates=0
- sigma=0.23: meshes=1,797, rows=7,188/7,188, failed=0, zero_probability=0, duplicates=0
- sigma=0.26: meshes=1,797, rows=7,188/7,188, failed=0, zero_probability=0, duplicates=0
- sigma=0.30: meshes=1,797, rows=7,188/7,188, failed=0, zero_probability=0, duplicates=0

## Overall residual metrics

- sigma=0.18: median=-0.149, mean=-0.165, MAE=0.381, RMSE=0.484, abs_p95=0.941
- sigma=0.20: median=-0.079, mean=-0.077, MAE=0.338, RMSE=0.422, abs_p95=0.801
- sigma=0.23: median=0.023, mean=0.044, MAE=0.300, RMSE=0.377, abs_p95=0.739
- sigma=0.26: median=0.112, mean=0.152, MAE=0.293, RMSE=0.381, abs_p95=0.769
- sigma=0.30: median=0.224, mean=0.281, MAE=0.341, RMSE=0.434, abs_p95=0.874

## Best sigma by probability level

- ALL: best RMSE sigma=0.23 (RMSE=0.377); best MAE sigma=0.26 (MAE=0.293)
- T50_P02_BV: best RMSE sigma=0.23 (RMSE=0.469); best MAE sigma=0.23 (MAE=0.377)
- T50_P05_BV: best RMSE sigma=0.23 (RMSE=0.397); best MAE sigma=0.26 (MAE=0.324)
- T50_P10_BV: best RMSE sigma=0.26 (RMSE=0.340); best MAE sigma=0.26 (MAE=0.281)
- T50_P39_BV: best RMSE sigma=0.30 (RMSE=0.218); best MAE sigma=0.26 (MAE=0.183)

## Interpretation

- The global RMSE minimum is sigma=0.23, while MAE selects sigma=0.26.
- The optimum is probability-level dependent, so a single fixed sigma should be described as an effective calibration parameter for the public-parameter prototype, not as a recovered proprietary J-SHIS internal parameter.
- The 0.20-0.26 range brackets the global minimum and supports using sigma=0.23 as the central production setting for the station-correction experiments, with sensitivity checks reported around it.
- This resolves the parameter-choice weakness by replacing an unexplained fixed value with a reproducible national stratified sensitivity audit.

Figure: `outputs/jshis_public_psha_sigma_sensitivity1797_rmse.png`

## Files

- `outputs/jshis_public_psha_sigma_sensitivity1797_metrics.csv`
- `outputs/jshis_public_psha_sigma_sensitivity1797_best.csv`
- `outputs/jshis_public_psha_sigma_sensitivity1797_integrity.csv`
- `outputs/jshis_public_psha_sigma_sensitivity1797_rows.csv`