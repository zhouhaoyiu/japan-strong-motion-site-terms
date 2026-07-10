# MF2013 magnitude, distance and event-station sensitivity

Morikawa and Fujiwara (2013) selected regression records with Mw >= 5.5, ground-surface sensors, at least five triggered stations per event and source distance below 200 km. This sensitivity applies the magnitude, distance and event-station-count conditions to the primary surface sample. It does not reconstruct the paper's additional magnitude-dependent PGA truncation.

- Records at each period: 53,517
- Earthquakes: 488
- Stations with at least 20 records: 772
- SA(3.0 s) station-field correlation with the primary field: 0.959
- SA(3.0 s) station-field 95th-percentile absolute difference: 0.083 log10 units
- SA(3.0 s) restricted event-holdout correlation: 0.895
- SA(3.0 s) restricted spatial RMSE gain: 9.9%
- SA(3.0 s) primary/restricted OOF-prediction correlation: 0.770
- SA(3.0 s) restricted multiplier 5th--95th percentiles: 0.649--1.445
- SA(3.0 s) restricted surface and adjusted medians: 0.089 and 0.088 g
- On the same 772 stations, the primary adjusted median is 0.091 g and the restricted adjusted median is 0.088 g, from a common surface median of 0.089 g

The restricted calculation tests selected magnitude, distance and event-station conditions from the original regression design. The primary analysis retains the broader public-flatfile domain because it evaluates the implementation used with the national response-spectrum product. Results outside the implemented screens are interpreted through this sensitivity.
