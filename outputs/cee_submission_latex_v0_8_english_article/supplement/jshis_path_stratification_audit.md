# Path and source stratification audit

## Design

- The primary RotD100 residuals are re-estimated within four hypocentral-bearing sectors, four shortest-fault-distance bins and three J-SHIS source classes.
- Each stratum is separated into connected station-event graph components and solved by sparse least squares. This preserves the primary additive event-station model without imposing offsets between disconnected components.
- A station requires at least 10 records within a stratum. A connected component requires at least 10 supported stations. Each retained component is aligned to the full-sample field by its record-weighted mean difference before comparison.
- Sampling reliability is evaluated with repeated disjoint event halves at 1, 2 and 3 s. Each station requires at least 5 records per half. We report both the reference-aligned correlation and a correlation after removing an independent weighted offset from every intersecting component pair.

## Results at 1-3 s

- Hypocentral-bearing sectors: weighted station-term correlations span 0.280--0.943; the largest weighted 95th-percentile absolute difference is 0.474 log10 units, with at least 448 paired stations.
- Shortest-fault-distance bins: weighted station-term correlations span 0.496--0.987; the largest weighted 95th-percentile absolute difference is 0.427 log10 units, with at least 138 paired stations.
- Source classes: weighted station-term correlations span 0.924--0.984; the largest weighted 95th-percentile absolute difference is 0.151 log10 units, with at least 960 paired stations.

## Split-half sampling reliability

- Hypocentral-bearing sectors: mean within-component split-half correlations span 0.764--0.969; each split retains at least 300 paired stations.
- Shortest-fault-distance bins: mean within-component split-half correlations span 0.631--0.968; each split retains at least 77 paired stations.
- Source classes: mean within-component split-half correlations span 0.902--0.981; each split retains at least 860 paired stations.
- Largest absolute change after removing component-pair offsets: 0.134.
- Maximum relative normal-equation residual across all component solves: 1.899e-09.
- Component solves outside accepted LSMR stop codes: 0.

## Interpretation boundary

The bearing is the initial great-circle direction from the event horizontal coordinates to the station; it is a reproducible directional proxy rather than a finite-fault ray path. Component-specific alignment removes offsets that are not identifiable in disconnected station-event graphs. Reported split-half correlations are invariant to those offsets because each intersecting component pair is centred independently. Distance and source-class subsets contain fewer records per station, so their differences combine sampling uncertainty with any unresolved path dependence. This audit tests whether the station field is dominated by one path sector or source class. It does not identify a path-specific nonergodic term or fully separate site and path effects.
