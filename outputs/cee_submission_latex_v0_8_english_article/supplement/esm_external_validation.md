# ESM external station-field validation

## Data boundary

- Source: Engineering Strong-Motion Database (ESM) v3 flatfile web service.
- Local snapshot SHA-256: `084ef4c7519fdeb26a05db9ad980535cd4c2207b2358d6b668cd5dafb4f0c724`.
- Snapshot rows: 41,841; events through 2026-07-07 18:11:32.
- The primary sample retains European-Mediterranean shallow events with 4.5 <= Mw <= 7.6 from 2015 onward, after publication of the fixed Bindi et al. (2014) backbone.
- As-recorded horizontal geometric-mean spectra match the primary backbone component. Kotha et al. (2020), its Mw <= 7.4 range and RotD50 are retained as a smaller response-component and backbone sensitivity.
- Event-station graph components are solved separately. Independent offsets are removed within intersecting train-test component pairs before pooled comparison.

## Primary post-2015 replication

- Supported field: 265 stations with at least 15 records across 3 periods.
- SA(1 s): mean held-event correlation 0.929 (station-bootstrap 95% CI 0.899 to 0.950); RMSE gain 66.8% (60.5% to 71.6%).
- SA(2 s): mean held-event correlation 0.939 (station-bootstrap 95% CI 0.915 to 0.957); RMSE gain 69.5% (64.7% to 73.6%).
- SA(3 s): mean held-event correlation 0.928 (station-bootstrap 95% CI 0.896 to 0.951); RMSE gain 66.4% (60.9% to 70.9%).

## Sensitivities

- Kotha/RotD50 SA(1 s): correlation 0.860; RMSE gain 48.8%.
- Kotha/RotD50 SA(2 s): correlation 0.916; RMSE gain 58.1%.
- Kotha/RotD50 SA(3 s): correlation 0.919; RMSE gain 58.5%.

## Cross-region transfer boundary

- A frozen Japan K-NET model using only log10(Vs30) and elevation gives SA(3.0 s) correlation 0.047 and RMSE gain -2.8% in Europe.
- The adverse transfer is retained. Japan-specific basin-depth and volcanic-front variables are unavailable in ESM, so this test evaluates only the two harmonised variables and does not invalidate the Japan regional model.
- The external result supports repeatability of long-period station fields across regions. It does not support a globally transferable station predictor.

## ESM local spatial prediction

- K-means requested 5 geographic blocks. The fixed minimum of 30 stations per block selected 3 blocks without using response values.
- Three feature sets were fixed before evaluation: harmonised VS30 and elevation; all available ESM site metadata; and the same metadata plus coordinates. Hyperparameters match the Japan analysis and are not tuned on ESM outcomes.
- The strict design estimates the training field from four event groups, predicts spatially held stations, and scores predictions against the fifth event group. Training targets are centred only on spatial-training stations; out-of-fold predictions are centred without response labels.
- full_field_spatial_blocks, common_site_hgb, SA(3.0 s): correlation 0.074 (95% CI -0.043 to 0.189); RMSE gain -3.6% (-9.8% to 1.3%).
- full_field_spatial_blocks, esm_metadata_hgb, SA(3.0 s): correlation 0.052 (95% CI -0.061 to 0.163); RMSE gain -5.5% (-12.8% to -0.0%).
- full_field_spatial_blocks, esm_metadata_location_hgb, SA(3.0 s): correlation 0.017 (95% CI -0.090 to 0.131); RMSE gain -7.6% (-13.1% to -2.2%).
- held_event_field_spatial_blocks, common_site_hgb, SA(3.0 s): correlation 0.021 (95% CI -0.083 to 0.132); RMSE gain -3.7% (-8.6% to 0.5%).
- held_event_field_spatial_blocks, esm_metadata_hgb, SA(3.0 s): correlation 0.001 (95% CI -0.095 to 0.097); RMSE gain -6.4% (-11.9% to -2.3%).
- held_event_field_spatial_blocks, esm_metadata_location_hgb, SA(3.0 s): correlation 0.006 (95% CI -0.094 to 0.102); RMSE gain -5.4% (-10.2% to -1.6%).
- All spatial results, including adverse folds and intervals, are retained. They test geographic extrapolation of available public ESM variables and do not test unavailable basin-depth variables.

## Selection table

```text
analysis	stage	n_records	n_events	n_stations
bindi2014_rhyp_mw4.5-7.6_2015-01-01	raw_query	41841	4194	3197
bindi2014_rhyp_mw4.5-7.6_2015-01-01	physical_and_response_filters	14463	853	959
bindi2014_rhyp_mw4.5-7.6_2015-01-01	event_station_deduplication	14463	853	959
bindi2014_rhyp_mw4.5-7.6_2015-01-01	iterative_event_station_support	13430	634	539
bindi2014_rhyp_mw5-7.6_2015-01-01	raw_query	41841	4194	3197
bindi2014_rhyp_mw5-7.6_2015-01-01	physical_and_response_filters	6041	283	776
bindi2014_rhyp_mw5-7.6_2015-01-01	event_station_deduplication	6041	283	776
bindi2014_rhyp_mw5-7.6_2015-01-01	iterative_event_station_support	5411	213	371
kotha2020_rjb_mw4.5-7.4_all	raw_query	41841	4194	3197
kotha2020_rjb_mw4.5-7.4_all	physical_and_response_filters	2075	139	490
kotha2020_rjb_mw4.5-7.4_all	event_station_deduplication	2075	139	490
kotha2020_rjb_mw4.5-7.4_all	iterative_event_station_support	1616	72	193
```
