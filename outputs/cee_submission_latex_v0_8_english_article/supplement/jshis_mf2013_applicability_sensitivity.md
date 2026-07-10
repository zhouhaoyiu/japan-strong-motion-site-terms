# MF2013 regression-domain sensitivity

Morikawa and Fujiwara (2013) selected regression records with Mw >= 5.5, ground-surface sensors, at least five triggered stations per event and source distance below 200 km. They also removed records beyond the distance at which the Kanno et al. (2006) median PGA fell below 10 cm/s/s. This sensitivity applies all of those conditions to the primary surface sample. The Kanno base equations use separate shallow and deep coefficients at a focal-depth boundary of 30 km; their station-amplification term is excluded from this magnitude-distance boundary.

- Records before and after the Kanno-PGA screen: 53,737 and 36,037
- Records at each period: 35,857
- Earthquakes: 411
- Stations with at least 20 records: 578
- SA(3.0 s) station-field correlation with the primary field: 0.959
- SA(3.0 s) station-field 95th-percentile absolute difference: 0.082 log10 units
- SA(3.0 s) restricted event-holdout correlation: 0.762
- SA(3.0 s) restricted spatial RMSE gain: 4.1%
- SA(3.0 s) primary/restricted OOF-prediction correlation: 0.648
- SA(3.0 s) restricted multiplier 5th--95th percentiles: 0.660--1.503
- SA(3.0 s) restricted surface and adjusted medians: 0.094 and 0.096 g
- On the same 578 stations, the primary adjusted median is 0.098 g and the restricted adjusted median is 0.096 g, from a common surface median of 0.094 g

The restricted calculation tests the published regression-domain conditions using the current public flatfile. It does not recreate the historical waveform database or the original regression weights. The primary analysis retains the broader public-flatfile domain because it evaluates the implementation used with the national response-spectrum product. Results outside the regression domain are interpreted through this sensitivity.
