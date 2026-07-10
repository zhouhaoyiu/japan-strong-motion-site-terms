# Network, geography, and influential-event stress tests

## Validation design

- Cross-network transfer trains the station model on one network and evaluates it on the other without using target-network station terms.
- Fixed macroregion tests leave out one of five coordinate-defined regions. These labels are reproducible geographic partitions, not administrative regions.
- The common-feature sensitivity excludes VS20, VS30 and network identity because their completeness differs between K-NET and KiK-net; it retains VS10, AVS30, sedimentary depths, elevation, sensor depth and geographic variables.
- Influential-event tests remove each of the 10 events with the largest SA(3.0 s) record counts and then remove all 10 jointly.
- Event-deletion comparisons align the perturbed and full station-term references by their record-weighted mean difference before computing changes.

## Results

- SA(3.0 s) cross-network weighted correlations span 0.417--0.619.
- SA(3.0 s) cross-network RMSE gains span 4.1%--21.2%.
- With common-coverage features, SA(3.0 s) cross-network correlations span 0.469--0.607, and RMSE gains span 8.7%--20.5%.
- SA(3.0 s) leave-macroregion-out RMSE gains span 5.4%--28.9%.
- Common-feature leave-macroregion-out SA(3.0 s) RMSE gains span -5.2%--27.5%.
- Individual high-count event deletion retains SA(3.0 s) station-term correlations of 0.9999--1.0000.
- Joint top-10 deletion retains an SA(3.0 s) correlation of 0.9994, with a weighted 95th-percentile absolute change of 0.012 log10 units.
- Joint top-10 deletion correlations across all periods span 0.9989--0.9998.

## Interpretation boundary

The transfer tests reuse station terms estimated from the same observational archive; they are model-generalization tests rather than independent-network replications. Event deletion tests sensitivity to high-count events and do not replace a fully independent earthquake catalogue.
