# Nearest-station official borehole-profile audit for the 10,467-mesh run

## Purpose

This audit adds an independent station-level physical layer using public NIED Kyoshin station-information profiles for the nearest K-NET/KiK-net station matched to each sampled mesh. It extracts upper-20 m borehole VS proxies and compares residual-tail, northeast-Hokkaido, and D1400 period-window groups. The audit is a nearest-station proxy check and does not replace three-dimensional basin-response modeling.

## Coverage

- Unique nearest stations requested: 1,543.
- Meshes linked to a parsable official station profile: 10,321 of 10,467 (98.6%).
- Residual-tail profile coverage: 1,739 of 1,743 (99.8%).
- Northeast-Hokkaido profile coverage: 2,553 of 2,553 (100.0%).

## Main findings

- Tail meshes with official profiles have a median borehole VS20 proxy of 284 m/s, compared with 316 m/s for non-tail meshes.
- Tail meshes have a median upper-20 m thickness with VS <=300 m/s of 6.0 m, compared with 5.0 m for non-tail meshes.
- Northeast-Hokkaido meshes with official profiles have a median borehole VS20 proxy of 253 m/s, compared with 327 m/s outside the band.
- Meshes in the D1400-derived 1--3 s window have a median borehole VS20 proxy of 287 m/s, compared with 322 m/s outside the window.

## Interpretation

The borehole-profile subset supports a shallow-velocity component in the residual-tail structure. Coverage is high overall, including complete coverage for the northeast-Hokkaido subset, but the audit remains a nearest-station proxy rather than a three-dimensional basin-response model. It should be used as independent consistency evidence, not as a national full-coverage geologic mechanism test.

## Output files

- `outputs/jshis_public_psha_10467_nearest_station_borehole_profiles.csv`
- `outputs/jshis_public_psha_10467_nearest_station_borehole_mesh_link.csv`
- `outputs/jshis_public_psha_10467_nearest_station_borehole_summary.csv`
- `outputs/cee_submission_latex_v0_8_english_article/figures/figure14_nearest_station_borehole_audit.png`
