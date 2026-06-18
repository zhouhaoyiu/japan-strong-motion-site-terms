# Basin-period proxy audit for the 10,467-mesh residual-tail structure

## Scope

This audit adds an independent physical-proxy layer to the public-parameter J-SHIS PBV residual analysis. It uses the 10,467-mesh nearest-station table, official D1400/Dbase/AVS30/VS20 site fields, and MF2013 site-corrected station residuals. D1400 period proxies are computed as 4D1400/700 and 4D1400/1000. These proxies indicate whether the nearby basin-depth scale is compatible with the 1--3 s band examined in the main manuscript.

## Main findings

- The highest basin long-period proxy quartile has a residual-tail fraction of 17.1%, compared with 16.5% in the lower quartiles; the tail risk ratio is 1.04.
- The northeast-Hokkaido band has a residual-tail fraction of 45.5%, compared with 7.3% outside the band; the tail risk ratio is 6.20.
- Meshes whose D1400 proxy lies in the 1--3 s window have a residual-tail fraction of 21.2%, compared with 14.3% outside that window.
- The median 4D1400/700 proxy is 1.21 s in northeast Hokkaido and 0.43 s outside that band.
- The same proxy is 0.95 s for residual-tail meshes and 0.69 s for non-tail meshes.
- Nearby station SA(3.0 s)-minus-PGA site-residual contrasts remain modest: the median is -0.341 for residual-tail meshes and -0.221 for non-tail meshes. The stronger separation comes from basin-depth and shallow-velocity proxies.

## Strongest monotonic associations with maximum absolute residual

- official_bv_max: Spearman rho 0.424; tail median 83.016; non-tail median 58.724.
- n_candidate_sources: Spearman rho -0.255; tail median 3295148.000; non-tail median 4423188.500.
- site_residual_sa3_minus_sa03: Spearman rho -0.206; tail median -0.263; non-tail median -0.158.
- site_residual_sa3_minus_pga: Spearman rho -0.194; tail median -0.341; non-tail median -0.221.
- basin_long_period_index: Spearman rho 0.058; tail median -0.346; non-tail median -0.481.

## Interpretation

The residual-tail structure is consistent with a basin-period proxy: the strongest geographic tail occurs where nearby official site fields imply deeper velocity horizons and longer D1400 period proxies. Station residual period profiles are a secondary signal in this audit. The evidence supports physical consistency between long-period site variables and the public-parameter PBV residual pattern; source aggregation and three-dimensional basin structure remain necessary for full mechanism attribution at individual 250 m meshes.

## Output files

- `outputs/jshis_public_psha_10467_basin_period_proxy_link.csv`
- `outputs/jshis_public_psha_10467_basin_period_proxy_summary.csv`
- `outputs/jshis_public_psha_10467_basin_period_proxy_associations.csv`
- `outputs/jshis_public_psha_10467_basin_period_proxy_tail_risk.csv`
- `outputs/cee_submission_latex_v0_8_english_article/figures/figure12_basin_period_proxy.png`
