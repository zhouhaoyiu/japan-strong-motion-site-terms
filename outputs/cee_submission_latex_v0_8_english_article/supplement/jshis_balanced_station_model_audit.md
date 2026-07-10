# Equal-stratum station-field sensitivity

## Design

- Component-aligned station terms are averaged with equal weight across available hypocentral-bearing sectors, distance bins or source classes.
- A station must contribute to at least 2 strata. The public-variable model uses the same five-block clustering procedure and hyperparameters as the primary analysis. Target strata receive equal weight; validation errors retain station record-count weights.

## Results at 1-3 s

- Azimuth-balanced field: 1,220 stations at 3 s; RMSE gains span 3.2%--6.4% and correlations span 0.306--0.406; the minimum held-block RMSE gain is -66.4%.
- Distance-balanced field: 1,285 stations at 3 s; RMSE gains span 2.5%--7.3% and correlations span 0.347--0.446; the minimum held-block RMSE gain is -9.6%.
- Source-class-balanced field: 1,369 stations at 3 s; RMSE gains span 5.0%--9.6% and correlations span 0.422--0.509; the minimum held-block RMSE gain is -24.5%.

## Interpretation boundary

Equal-stratum averaging reduces dominance by the observed path distribution. It does not identify a pure site term, because station-specific path interactions can remain within each stratum and coverage is incomplete. These results are a sensitivity analysis of the persistent station-associated component.
