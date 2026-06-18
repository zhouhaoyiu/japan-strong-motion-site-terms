# GSI terrain geomorphic proxy audit

This audit samples Geospatial Information Authority of Japan (GSI) public DEM PNG elevation tiles at zoom 10 for the 10,467 sampled public-parameter meshes. It is independent of MF2013 site terms and J-SHIS station residuals. It supplies a geomorphic proxy check, not a three-dimensional velocity model.

## Inputs

- Mesh table: `outputs/cee_submission_latex_v0_8_english_article/supplement/jshis_sa3_continuous_correction_surface_mesh.csv`.
- GSI layer: `dem_png`.
- Tile cache: `work/external_data/gsi_dem_png_z10`.
- New tile downloads in this run: 0.

## Coverage

- Sampled meshes: 10,467.
- Valid GSI elevation samples: 9,956 (95.1%).
- Median GSI elevation: 168.4 m.
- Median 1 km relief: 114.0 m.

## Main contrasts

- Residual-tail elevation median: 148.1 m; non-tail median: 174.1 m.
- Northeast-Hokkaido elevation median: 153.5 m; outside-northeast-Hokkaido median: 174.2 m.
- Residual-tail lowland-flat fraction: 15.1%; non-tail fraction: 13.4%.
- Northeast-Hokkaido lowland-flat fraction: 13.7%; outside-northeast-Hokkaido fraction: 13.7%.
- Lowland-flat tail fraction: 18.4%; other-terrain tail fraction: 16.5%; risk ratio: 1.12.
- Tail versus non-tail elevation rank-biserial statistic: -0.061; p value: 8.82e-05.

## Boundary

The GSI terrain fields describe surface geomorphology. They can test whether residual-tail and northeast-Hokkaido meshes align with low-elevation or low-relief terrain expected for sedimentary settings. They do not measure shear-wave velocity, basin-edge geometry, or nonlinear site response.
