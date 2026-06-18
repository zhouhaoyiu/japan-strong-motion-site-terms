# J-SHIS national stratified 5322-mesh 12-process recomputation audit

- Shards: 12
- Sample meshes: 5,322
- Probability rows: 21,288
- Expected probability rows: 21,288
- Status counts: `{'ok': 21288}`
- Empty/failed rows: 0
- Zero-probability rows: 0
- Duplicate mesh/probability rows: 0
- Sum of tile compute time: 17804.54 s
- Approximate wall-clock compute time from slowest shard: 1495.49 s
- Effective wall-clock seconds per mesh: 0.281
- Mean tile compute seconds per mesh, summed over processes: 3.345
- Wall-clock mode: 12 parallel single-core processes
- Candidate sources median by tile: 4105640
- Candidate sources maximum tile mean: 4794337

## Residual summary at official PBV thresholds

- T50_P02_BV: median=0.006, p05=-0.540, p95=0.846, min=-1.199, max=1.113
- T50_P05_BV: median=0.007, p05=-0.566, p95=0.675, min=-0.918, max=0.905
- T50_P10_BV: median=-0.009, p05=-0.601, p95=0.555, min=-0.911, max=0.733
- T50_P39_BV: median=-0.019, p05=-0.624, p95=0.271, min=-0.781, max=0.334

## Candidate/source diagnostics

- n_candidate_sources: median=4215614, p05=2284633, p95=4867777, max=4898504
- n_sources_used: median=4099000, p05=2183576, p95=4708498, max=4737430
- n_geometry_points_used: median=4209914, p05=2275387, p95=4863315, max=4894693

## Parallel shard timing

- national_d1200_sigma0p23_stratified5322_arrays12_shard00: sites=444, compute=1450.10 s, max_tile_mean_candidates=4785536
- national_d1200_sigma0p23_stratified5322_arrays12_shard01: sites=444, compute=1495.49 s, max_tile_mean_candidates=4789463
- national_d1200_sigma0p23_stratified5322_arrays12_shard02: sites=444, compute=1490.13 s, max_tile_mean_candidates=4794337
- national_d1200_sigma0p23_stratified5322_arrays12_shard03: sites=444, compute=1493.77 s, max_tile_mean_candidates=4784177
- national_d1200_sigma0p23_stratified5322_arrays12_shard04: sites=444, compute=1484.58 s, max_tile_mean_candidates=4786199
- national_d1200_sigma0p23_stratified5322_arrays12_shard05: sites=444, compute=1472.34 s, max_tile_mean_candidates=4792684
- national_d1200_sigma0p23_stratified5322_arrays12_shard06: sites=443, compute=1474.64 s, max_tile_mean_candidates=4776782
- national_d1200_sigma0p23_stratified5322_arrays12_shard07: sites=443, compute=1486.66 s, max_tile_mean_candidates=4777718
- national_d1200_sigma0p23_stratified5322_arrays12_shard08: sites=443, compute=1492.69 s, max_tile_mean_candidates=4776299
- national_d1200_sigma0p23_stratified5322_arrays12_shard09: sites=443, compute=1495.09 s, max_tile_mean_candidates=4777153
- national_d1200_sigma0p23_stratified5322_arrays12_shard10: sites=443, compute=1488.56 s, max_tile_mean_candidates=4776646
- national_d1200_sigma0p23_stratified5322_arrays12_shard11: sites=443, compute=1480.49 s, max_tile_mean_candidates=4788591

This audit is a deterministic national stratified sample over geographic bins, hazard quantiles, and geographic extremes. It verifies national-scale tiled execution under 12-process parallelism. It is not the completed 5,989,018-mesh national run.
