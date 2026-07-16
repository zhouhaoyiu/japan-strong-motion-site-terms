# Source-conditioned station-field and response-spectrum sensitivity

## Locked design

- Source classes are the public flatfile labels: crustal, interplate and intraplate. Every class, period and spatial fold is retained.
- A station requires at least 20 records within a source class. Targets are record-weighted to zero within source class and period before fitting.
- The fixed primary site-plus-location model and 5 spatial blocks are used without tuning. Confidence intervals use 2,000 station-cluster bootstrap replicates.
- The crustal field is applied to the official active-shallow product. Interplate and intraplate fields are applied separately to the official subduction product as two sensitivities; they are not averaged because official source weights are unavailable.

## Spatial prediction

### Crustal

- 1 s: 1,084 stations; correlation 0.484 (0.438--0.530); RMSE gain 12.7% (8.2%--17.2%).
- 2 s: 1,084 stations; correlation 0.468 (0.414--0.517); RMSE gain 10.5% (5.8%--14.8%).
- 3 s: 1,084 stations; correlation 0.455 (0.403--0.507); RMSE gain 8.9% (3.9%--13.8%).
- 5 s: 1,084 stations; correlation 0.383 (0.329--0.437); RMSE gain 4.2% (-0.7%--8.9%).

### Interplate

- 1 s: 779 stations; correlation 0.384 (0.323--0.444); RMSE gain 4.5% (-1.1%--9.9%).
- 2 s: 779 stations; correlation 0.228 (0.153--0.301); RMSE gain -0.8% (-6.2%--4.5%).
- 3 s: 779 stations; correlation 0.251 (0.175--0.322); RMSE gain 2.8% (-1.7%--7.4%).
- 5 s: 779 stations; correlation 0.298 (0.228--0.360); RMSE gain 1.7% (-2.8%--5.9%).

### Intraplate

- 1 s: 856 stations; correlation 0.404 (0.345--0.461); RMSE gain 4.6% (-0.1%--9.1%).
- 2 s: 856 stations; correlation 0.396 (0.332--0.453); RMSE gain 7.4% (3.0%--11.6%).
- 3 s: 856 stations; correlation 0.516 (0.460--0.567); RMSE gain 14.3% (10.1%--18.3%).
- 5 s: 856 stations; correlation 0.548 (0.502--0.594); RMSE gain 14.8% (10.7%--18.6%).

## Response-spectrum propagation

- active_shallow_crustal: 1,084 matched stations; median official 0.033 g; median surface reference 0.034 g; median adjusted 0.033 g.
- subduction_interplate: 779 matched stations; median official 0.077 g; median surface reference 0.082 g; median adjusted 0.084 g.
- subduction_intraplate: 856 matched stations; median official 0.082 g; median surface reference 0.089 g; median adjusted 0.083 g.

## Interpretation boundary

These fields condition the empirical station residual on the observed source class. They retain unresolved path effects within each class and are sensitivities on published fixed-probability products. They do not replace source-specific terms inside a production PSHA integral.
