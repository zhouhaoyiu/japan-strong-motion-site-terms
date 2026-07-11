# Cross-network and held-event station-term validation

K-NET station terms are estimated only from training events. The fixed public-variable model is then evaluated against KiK-net station terms estimated only from disjoint held-out events. The reverse direction is retained as a transfer stress test.

- Event folds: 5; fixed seed: 20260711.
- SA(3.0 s) K-NET to KiK-net target stations: 677.
- SA(3.0 s) Pearson correlation: 0.598 (station bootstrap 95% CI 0.542 to 0.651).
- SA(3.0 s) RMSE gain: 20.8% (95% CI 13.3% to 26.7%).
- Reverse KiK-net to K-NET correlation and RMSE gain: 0.462 and 4.0%.

The validation separates both recording network and earthquake set. No target-network station term, target event, or target prediction is used to fit or recenter the frozen model.
