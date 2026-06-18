# KiK-net Surface/Downhole Amplification Proxy Audit

## Data inventory
- Paired surface/downhole miniSEED records: 288
- Sampling rates: [100.0]
- Median record length: 24,250 samples

## Horizontal amplification proxy
- Median surface/downhole horizontal PGA ratio: 1.111
- Interquartile range: 0.650 to 1.614
- Stations with ratio > 1: 161 / 288
- Stations with ratio > 2: 39 / 288
- Median horizontal RMS ratio: 1.146
- Spearman correlation with downhole-motion level: -0.587 (single-event diagnostic only)

## Frequency-dependent ratios
- 0.2-0.5 Hz median spectral ratio: 2.382
- 0.5-1 Hz median spectral ratio: 1.857
- 1-2 Hz median spectral ratio: 1.152
- 2-5 Hz median spectral ratio: 0.667
- 5-10 Hz median spectral ratio: 0.464
- 10-20 Hz median spectral ratio: 0.694

## Highest horizontal PGA ratios
- TYMH04: PGA ratio=2.76, horizontal_rms_ratio=2.30, SR 1_2Hz=1.80, SR 2_5Hz=2.22, SR 5_10Hz=1.55
- AICH15: PGA ratio=2.62, horizontal_rms_ratio=2.60, SR 1_2Hz=2.22, SR 2_5Hz=1.10, SR 5_10Hz=0.35
- AICH19: PGA ratio=2.58, horizontal_rms_ratio=2.62, SR 1_2Hz=2.31, SR 2_5Hz=1.58, SR 5_10Hz=1.75
- KKWH03: PGA ratio=2.55, horizontal_rms_ratio=2.44, SR 1_2Hz=1.86, SR 2_5Hz=1.06, SR 5_10Hz=0.78
- SZOH53: PGA ratio=2.54, horizontal_rms_ratio=2.57, SR 1_2Hz=2.18, SR 2_5Hz=1.12, SR 5_10Hz=1.14
- NGNH54: PGA ratio=2.54, horizontal_rms_ratio=2.39, SR 1_2Hz=2.01, SR 2_5Hz=0.81, SR 5_10Hz=0.60
- SOYH08: PGA ratio=2.51, horizontal_rms_ratio=2.50, SR 1_2Hz=1.68, SR 2_5Hz=0.36, SR 5_10Hz=0.36
- SZOH30: PGA ratio=2.50, horizontal_rms_ratio=2.55, SR 1_2Hz=2.33, SR 2_5Hz=1.33, SR 5_10Hz=0.38
- GIFH16: PGA ratio=2.47, horizontal_rms_ratio=2.47, SR 1_2Hz=2.02, SR 2_5Hz=1.13, SR 5_10Hz=1.52
- NIGH15: PGA ratio=2.45, horizontal_rms_ratio=2.06, SR 1_2Hz=1.46, SR 2_5Hz=0.78, SR 5_10Hz=0.85
- NGNH29: PGA ratio=2.39, horizontal_rms_ratio=1.76, SR 1_2Hz=0.78, SR 2_5Hz=0.55, SR 5_10Hz=0.41
- GIFH19: PGA ratio=2.38, horizontal_rms_ratio=2.50, SR 1_2Hz=2.13, SR 2_5Hz=1.09, SR 5_10Hz=1.06

## Interpretation guardrails
- These are surface/downhole amplification proxies from one event and count-like waveforms after detrending.
- They support candidate-site discovery and frequency-band targeting, not a final nonlinear site-response claim.
- The next publication-grade step is to merge full KiK-net multi-event data with station metadata such as borehole depth and Vs profiles.
