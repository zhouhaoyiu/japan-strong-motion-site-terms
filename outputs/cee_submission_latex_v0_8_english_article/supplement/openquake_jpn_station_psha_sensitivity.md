# OpenQuake Japan PSHA station-term sensitivity

This is an independent OpenQuake/GEM Japan-model sensitivity run for five representative stations. It is not an official J-SHIS hazard-map reproduction.

The OpenQuake calculation uses the packaged Japan mosaic source-model logic tree and GMPE logic tree. The station term is applied after the baseline calculation as a deterministic SA(3.0 s) spectral shift. For station multiplier f, lambda_c(S)=lambda_0(S/f).

## 50-year 10% SA(3.0 s) summary

| site_code | example_type | station_multiplier | baseline_g | corrected_g | delta_pct | official_jshis_50y10pct_sa3_g | oq_to_jshis_official_sa3_ratio_50y10pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| KYTH07 | maximum_negative | 0.1625 | 0.1623 | 0.0264 | -83.7453 | 0.1444 | 1.1238760080138823 |
| NIG024 | maximum_positive | 1.5498 | 0.0461 | 0.0715 | 54.9784 | 0.051 | 0.9050251445387962 |
| KOCH10 | q05 | 0.4136 | 0.1346 | 0.0557 | -58.6373 | 0.1497 | 0.898811411413917 |
| FKS031 | q50 | 0.6191 | 0.0183 | 0.0113 | -38.0945 | 0.0728 | 0.25073513879085496 |
| MYZH15 | q95 | 0.8619 | 0.1096 | 0.0944 | -13.8122 | 0.0753 | 1.455310163128457 |

## Interpretation

- The baseline OpenQuake PSHA produces full hazard curves and UHS ordinates for the selected stations.
- Applying the deterministic station multiplier changes only the SA(3.0 s) ordinate and leaves other periods fixed.
- The calculation is useful as an independent PSHA sensitivity check, while the manuscript's main official-data result remains the J-SHIS response-spectrum sensitivity audit.
