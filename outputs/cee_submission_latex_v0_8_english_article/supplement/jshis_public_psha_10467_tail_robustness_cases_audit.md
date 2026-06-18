# Tail robustness and case-profile audit for the 10,467-mesh run

## Purpose

This audit tests whether the residual-tail structure survives practical confounder checks. It stratifies by nearest-station distance, official PBV level, candidate-source load, and low-latitude island membership; fits a logistic tail model; and builds matched tail/non-tail mesh profiles for direct inspection.

## Main findings

- Within 50 km of a matched station, northeast-Hokkaido meshes have a tail fraction of 48.0% compared with 7.3% outside the band.
- The D1400-derived 1--3 s proxy window has tail fractions of 21.2% versus 14.3% outside the window across all meshes, and 21.7% versus 14.2% within the 50 km station-distance subset.
- Inside the middle official-PBV deciles q3--q8, the same 1--3 s proxy contrast is 26.6% versus 17.2%, reducing dependence on the highest hazard decile.
- Logistic adjustment for official PBV, candidate-source count, nearest-station distance, and station residual contrast gives an adjusted odds ratio of 26.20 for northeast Hokkaido and 0.93 for the 1--3 s D1400 proxy window.
- Across tail thresholds from 0.50 to 0.90 absolute log10 units, the northeast-Hokkaido tail-risk ratio ranges from 2.63 to 29.72.
- The case-profile table lists 9 tail cases and 9 matched non-tail controls, with probability-level residuals, official PBV values, source-candidate counts, station distance, site proxies, and period proxies.

## Interpretation

The northeast-Hokkaido tail remains strong inside the 50 km station-distance subset. The D1400-derived 1--3 s proxy contrast remains visible after excluding the highest official PBV decile and after restricting to near-station meshes. In the adjusted model, northeast Hokkaido and official PBV dominate the tail flag; the D1400 window does not add a significant independent term after those controls. The evidence supports a regional source-site-basin interaction and does not reduce the tail to sampling distance, a single hazard-amplitude bin, or one basin-depth threshold.

## Output files

- `outputs/jshis_public_psha_10467_tail_robustness_distance_bins.csv`
- `outputs/jshis_public_psha_10467_tail_robustness_stratified.csv`
- `outputs/jshis_public_psha_10467_tail_threshold_sensitivity.csv`
- `outputs/jshis_public_psha_10467_tail_robustness_logistic.csv`
- `outputs/jshis_public_psha_10467_tail_case_profiles.csv`
- `outputs/jshis_public_psha_10467_tail_case_probability_rows.csv`
- `outputs/cee_submission_latex_v0_8_english_article/figures/figure13_tail_robustness_cases.png`
