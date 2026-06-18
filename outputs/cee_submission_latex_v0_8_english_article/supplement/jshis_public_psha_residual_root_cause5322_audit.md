# Residual root-cause audit for the 5,322-mesh national stratified recomputation

This audit targets the main review risk left after the national recomputation and sigma sensitivity checks: whether the residual tails are random numerical noise or structured by geography, official hazard level, probability level, or source-candidate load.

## Core residual state

- Meshes: 5,322
- Probability rows: 21,288
- Overall median residual: -0.005
- Overall mean residual: 0.028
- Overall MAE: 0.306
- Overall RMSE: 0.381
- Rows with positive residual, where the public prototype gives higher probability at the official PBV: 49.5%
- Rows with negative residual, where the public prototype gives lower probability at the official PBV: 50.5%
- Rows with |log10 residual| >= 0.70: 6.2%
- Sites with at least one probability level at |log10 residual| >= 0.70: 17.3%

## Main structured effects

- Largest probability-level RMSE: T50_P02_BV with RMSE=0.467.
- Largest coarse geographic p95 site tail: northeast_hokkaido_band with p95 max |residual|=1.022 and 749 tail sites.
- Largest hazard-decile p95 site tail: q8 with p95 max |residual|=1.061 and 136 tail sites.
- Strongest monotonic associations with absolute residual:
  - official_bv_cm_s: Spearman rho=0.499
  - official_bv_log10: Spearman rho=0.499
  - lon: Spearman rho=0.318

## Root-cause interpretation

- The residual tails are structured enough to report as model-discrepancy diagnostics, not numerical failures. All rows are valid and no zero-probability rows occur in the 5,322-mesh audit.
- The probability-level dependence means that one scalar accuracy number is insufficient. The paper should report the long-return and short-return levels separately.
- The coarse geographic and hazard-decile diagnostics identify where the public-parameter prototype departs most from official maps. This supports a transparent limitation statement and gives a concrete target for future full production reproduction.
- These findings should be used to strengthen the CEE submission by showing where the method works, where it deviates, and how the deviations were audited.

## Review relevance

- This audit increases review robustness because it answers the likely reviewer question: whether the residual tails have been inspected and bounded.
- It is more valuable for the scientific argument than another small sample-size increase alone because it shows where the public-parameter calculation agrees with official products and where structured discrepancies remain.
- Combined with the 5,322-mesh correctness audit and five-point sigma sensitivity audit, this supports a transparent limitation statement and gives concrete targets for future production-grade reproduction.

## Output files

- `outputs/jshis_public_psha_residual_root_cause5322_by_probability.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_by_coarse_geo.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_by_selection_group.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_by_official_bv_quantile.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_by_candidate_quantile.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_site_by_coarse_geo.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_site_by_hazard_quantile.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_spearman.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_tail_rows_abs_ge_0p70.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_top_sites.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_bin_table.csv`
- `outputs/jshis_public_psha_residual_root_cause5322_summary.png`
