# J-SHIS VS400 site-amplification layer audit

This audit joins the 10,467 sampled public-parameter meshes to the public J-SHIS V4 surface-ground data file `Z-V4-JAPAN-AMP-VS400_M250.csv`. The joined fields are the engineering geomorphology code, gridded AVS30, and the peak-velocity amplification factor from engineering bedrock with Vs=400 m/s to the ground surface. The layer is independent of the MF2013 residual calculation and the hazard-curve residual rows, while still being a J-SHIS/NIED public site-amplification product.

## Inputs

- Mesh table: `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_gsi_terrain_geomorphic_mesh.csv`.
- J-SHIS site-amplification table: `work/jshis_full_psha_audit/extracted/zamp/Z-V4-JAPAN-AMP-VS400_M250/Z-V4-JAPAN-AMP-VS400_M250.csv`.
- J-SHIS V4 surface-ground documentation identifies IDAVS as average S-wave velocity in the upper 30 m and IDARV2 as the maximum velocity amplification factor from engineering bedrock to the surface.

## Coverage

- Sampled meshes: 10,467.
- Valid VS400 site-amplification joins: 10,467 (100.0%).
- Median gridded AVS30: 495.1 m/s.
- Median VS400 amplification factor: 0.8338.

## Main contrasts

- Residual-tail AVS30 median: 433.8 m/s; non-tail median: 510.4 m/s.
- Residual-tail VS400 amplification median: 0.9332; non-tail median: 0.8125.
- Northeast-Hokkaido AVS30 median: 431.1 m/s; outside-northeast-Hokkaido median: 510.4 m/s.
- Northeast-Hokkaido VS400 amplification median: 0.9382; outside-northeast-Hokkaido median: 0.8125.
- D1400-derived 1--3 s window AVS30 median: 429.1 m/s; outside-window median: 517.5 m/s.
- Lowest-quartile AVS30 meshes have a residual-tail fraction of 21.5%, compared with 14.9% outside that quartile; risk ratio: 1.44.
- Tail versus non-tail AVS30 rank-biserial statistic: -0.090; p value: 1.77e-09.
- Highest tail-rate geomorphology codes with n>=30: JCODE 9: n=438, tail=45.9%, AVS30=282.1 m/s; JCODE 16: n=254, tail=24.0%, AVS30=260.2 m/s; JCODE 6: n=233, tail=23.6%, AVS30=405.6 m/s.

## Boundary

This audit supplies gridded shallow-velocity and surface-amplification evidence. It is stronger than a terrain-only proxy, but it remains a surface-ground product and does not resolve three-dimensional basin-edge geometry, low-frequency amplification spectra, nonlinear response, or full non-ergodic PSHA propagation.
