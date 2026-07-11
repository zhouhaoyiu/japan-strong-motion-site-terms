# Chronological cross-network station-field validation

Events are ordered by JMA origin time. K-NET station terms are estimated from the earlier event block, and the frozen public-variable model is evaluated against KiK-net station terms estimated from the later block. Recording network, event set, and target stations are all external to model fitting.

- Prespecified early/late event splits: 70/30, 80/20, and 90/10.
- Primary 80/20 cutoff: training through 2016-11-22; testing from 2016-11-24.
- Primary SA(3.0 s) target: 547 KiK-net stations from 256 later events.
- Primary SA(3.0 s) Pearson correlation: 0.536 (station bootstrap 95% CI 0.469 to 0.600).
- Primary SA(3.0 s) RMSE gain: 19.4% (95% CI 13.6% to 24.9%).
- Across all eight periods in the primary split, correlations range from 0.460 to 0.635; RMSE gains range from 12.6% to 21.5%, and every 95% lower bound remains above zero.
- Across the three splits and 1-3 s periods, correlations range from 0.327 to 0.543; RMSE gains range from 12.6% to 19.4%.
- Reverse KiK-net-to-K-NET transfer at 3.0 s gives a correlation of 0.302 and an RMSE gain of -3.3% (95% CI -7.8% to 1.5%).

The target-network terms and predictions are not used for model fitting or recentering. Each chronological split and the adverse reverse-direction result are reported, so the result does not depend on selecting a favorable cutoff year or transfer direction.
