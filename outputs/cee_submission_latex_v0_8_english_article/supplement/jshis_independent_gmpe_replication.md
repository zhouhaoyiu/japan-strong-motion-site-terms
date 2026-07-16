# Independent GMPE replication

OpenQuake engine version: 3.25.1.

## Input mapping

- Zhao et al. (2006) active-crust, subduction-interface and subduction-slab variants are selected from the J-SHIS source class.
- The documented J-SHIS shortest fault distance is supplied as rupture distance.
- J-SHIS AVS30 is supplied as Vs30 for every station.
- Observed RotD50 spectral acceleration and Zhao predictions are compared in log10(cm/s2).

## Coverage

- 0.1 s: 222,664 records.
- 0.2 s: 222,664 records.
- 0.3 s: 222,664 records.
- 0.5 s: 222,664 records.
- 1 s: 222,664 records.
- 2 s: 222,664 records.
- 3 s: 222,664 records.
- 5 s: 222,664 records.

## SA(3.0 s) result

- MF2013-Zhao station-term Pearson correlation: 0.770.
- MF2013-Zhao station-term Spearman correlation: 0.775.
- Event-holdout station correlation: 0.950.
- Spatial-block RMSE reduction: 37.7%.

## Boundary

This analysis changes the regional ground-motion model while retaining the same observations. It measures model dependence of the station terms. Independent-network evidence is evaluated separately, and the Zhao calculation uses the available flatfile inputs.

## Convergence

- Maximum final fixed-effect change: 9.980e-11 log10 units.
- Maximum iteration count: 1105.
