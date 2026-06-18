# J-SHIS national stratified 10467-mesh 8-process recomputation audit

- Shards: 8
- Sample meshes: 10,467
- Probability rows: 41,868
- Expected probability rows: 41,868
- Status counts: `{'ok': 41868}`
- Empty/failed rows: 0
- Zero-probability rows: 0
- Duplicate mesh/probability rows: 0
- Sum of tile compute time: 59588.03 s
- Approximate wall-clock compute time from slowest shard: 7465.56 s
- Effective wall-clock seconds per mesh: 0.713
- Mean tile compute seconds per mesh, summed over processes: 5.693
- Wall-clock mode: 8 parallel single-core array-backend processes
- Candidate sources median by tile: 4183976
- Candidate sources maximum tile mean: 4848429

## Residual summary at official PBV thresholds

- T50_P02_BV: median=0.011, p05=-0.521, p95=0.838, min=-1.199, max=1.117
- T50_P05_BV: median=0.011, p05=-0.568, p95=0.669, min=-0.918, max=0.909
- T50_P10_BV: median=-0.012, p05=-0.609, p95=0.545, min=-0.911, max=0.736
- T50_P39_BV: median=-0.030, p05=-0.626, p95=0.268, min=-0.808, max=0.335

This audit is a deterministic national stratified sample over geographic bins,
hazard quantiles, and geographic extremes. It verifies national-scale tiled
execution at a larger sample size than the preceding audit. It is not the
completed 5,989,018-mesh national production run.
