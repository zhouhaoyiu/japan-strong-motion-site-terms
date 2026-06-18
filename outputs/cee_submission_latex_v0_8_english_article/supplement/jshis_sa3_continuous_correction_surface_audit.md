# SA(3.0 s) continuous correction-surface audit

This audit estimates station-level SA(3.0 s) correction terms on the 10,467-mesh public-parameter audit grid. The primary surface uses public site and location predictors with a gradient-boosted station model. A nearest-station inverse-distance-weighted surface is retained as a spatial-smoothing baseline. The audit is a sampled-grid continuity check, not a production non-ergodic PSHA implementation.

## Inputs

- Station correction table: `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_station_hazard_impact_by_station.csv`.
- Site schema for station coordinates: `outputs/jshis_site_schema_v2024_sub1.csv`.
- Sampled mesh table: `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_public_psha_10467_nearest_station_borehole_mesh_link.csv`.
- Stations with usable correction and coordinates: 2,261.
- Sampled meshes with usable coordinates: 10,467.

## Method

The primary surface uses longitude, latitude, station elevation, sensor depth, volcanic-front distances, shallow-velocity fields, and basin-depth fields. On the sampled mesh grid, non-coordinate fields are nearest-station proxies from public K-NET/KiK-net metadata. The IDW baseline uses the nearest 24 stations in an approximate kilometre coordinate system. Validation leaves out the stored spatial station blocks and predicts held-out station corrections from the remaining blocks. Metrics use station record counts as weights.

## Station validation

- Weighted RMSE: 0.0678 log10 units.
- Weighted MAE: 0.0515 log10 units.
- RMSE reduction versus zero additional correction: 70.7%.
- RMSE reduction versus a single national mean correction: 32.5%.
- IDW baseline RMSE reduction versus zero additional correction: 53.0%.
- IDW baseline RMSE reduction versus a single national mean correction: -8.3%.
- Median nearest training-station distance in held-out validation: 144.6 km.

## Sampled-grid surface

- Primary all-mesh multiplier 5th/50th/95th percentiles: 0.513, 0.658, 0.912.
- IDW all-mesh multiplier 5th/50th/95th percentiles: 0.513, 0.617, 0.747.
- Residual-tail median multiplier: 0.651; non-tail median multiplier: 0.659.
- Northeast-Hokkaido median multiplier: 0.638; outside-northeast-Hokkaido median multiplier: 0.662.
- Median nearest correction-station distance on sampled meshes: 7.9 km.

## Boundary

The output supports the path from station multipliers to a spatially continuous site-term prototype over the sampled audit grid. The mesh values use nearest-station site proxies, so they should be read as an auditable bridge between station terms and gridded hazard products. They do not replace a full PSHA recomputation with gridded site fields, source/path uncertainty propagation, and all-period logic-tree integration.
