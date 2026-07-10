# Spatial-block station-term uncertainty propagation

## Design

- Target interval: weighted empirical 90% prediction interval (0.05--0.95 error quantiles).
- Calibration: for each held spatial fold, signed out-of-fold errors are estimated only from the other folds.
- Weighting: station errors are weighted by the number of contributing strong-motion records.
- Propagation: lower, point, and upper station terms multiply the same matched J-SHIS surface-reference spectrum.

## Results

- Station-period intervals: 13,024 rows.
- Propagated station-period-probability values: 52,096 rows.
- SA(3.0 s) station coverage: 89.0%.
- SA(3.0 s) record-weighted coverage: 88.8%.
- SA(3.0 s) median upper/lower multiplier span: 3.00.
- Coverage range across eight periods: 86.7%--89.9%.

## Interpretation boundary

These intervals describe prediction error for the station-term model under the existing spatial-block validation design. They do not include uncertainty in earthquake occurrence, source characterization, the MF2013 backbone, J-SHIS hazard calculations, or the surface-reference conversion. The propagated bounds are conditional sensitivity coordinates, not complete non-ergodic PSHA confidence intervals.
