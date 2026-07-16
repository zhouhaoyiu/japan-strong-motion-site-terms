# KiK-net surface/borehole spectral validation

Surface and borehole RotD100 spectra are paired by earthquake and KiK-net site code. Transfer ratios are estimated from disjoint training events and used to predict station terms fitted to held-out surface records.

- Paired records at each period: 102,428
- Paired sites and earthquakes: 699 and 1,410
- SA(3.0 s) station-term/transfer Pearson and Spearman correlations: 0.315 and 0.415
- SA(3.0 s) Pearson correlation after log-borehole-depth adjustment: 0.327
- SA(3.0 s) transfer train/test correlation across event folds: 0.988
- SA(3.0 s) cross-event station-term correlation: 0.290
- SA(3.0 s) cross-event RMSE gain over a zero station term: 2.7%
- Two-way-centred eight-period spectral-shape correlation: 0.677 (cluster-bootstrap 95% CI 0.648 to 0.704).
- Median within-station spectral-shape correlation and positive fraction: 0.769 and 91.5%.
- Mean spectral-shape correlation with transfer ratios and station terms estimated from disjoint events: 0.635.

The paired-sensor result provides an independent physical correlate of the event-adjusted surface station terms. It validates their station-level period dependence. A nonlinear site-response model would require input-motion-dependent calibration beyond this comparison.
