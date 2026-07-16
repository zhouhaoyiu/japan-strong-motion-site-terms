# Sung et al. (2025) Kanto SA(5 s) external comparison

## Locked design

- The external field is the unmodified valid portion of the publisher-supplied Kanto T=5 s adjusted-site-term CSV. No geographic response-value screen is applied.
- Each station is assigned its nearest valid simulation cell only when the distance is at most 1.5 times the median grid spacing (1.201 km).
- Sung et al. report natural-log units. Japanese log10 station terms and out-of-fold predictions are multiplied by ln(10). Each matched field pair is centred independently because the two residual coordinates have arbitrary offsets.
- The source-matched crustal observed field is the primary physical-comparability test. Crustal out-of-fold predictions, the all-source field and the physical-only model are retained as diagnostics.
- Spatial-block intervals use 2,000 replicates and the published Kanto VCM correlation length of 31.7 km. Stations receive equal weight.
- Matching and inference settings were fixed before calculating cross-field associations.

## Results

The source-matched comparison contains 364 stations in 88 spatial blocks. Its Pearson correlation is 0.304 (0.163--0.434), and its Spearman correlation is 0.372 (0.227--0.494). Centred signs agree at 62.9% of stations (56.5%--69.5%).

| Field | Stations | Blocks | Pearson | 95% interval | Spearman | 95% interval | Sign agreement |
|---|---:|---:|---:|---:|---:|---:|---:|
| All-source observed | 389 | 88 | 0.122 | -0.009--0.258 | 0.173 | 0.021--0.332 | 55.8% |
| All-source spatial OOF | 389 | 88 | 0.349 | 0.243--0.446 | 0.407 | 0.280--0.513 | 63.5% |
| All-source physical OOF | 389 | 88 | 0.298 | 0.196--0.395 | 0.326 | 0.194--0.442 | 57.8% |
| Crustal observed | 364 | 88 | 0.304 | 0.163--0.434 | 0.372 | 0.227--0.494 | 62.9% |
| Crustal spatial OOF | 364 | 88 | 0.135 | 0.038--0.231 | 0.141 | 0.005--0.276 | 54.9% |
| Crustal physical OOF | 364 | 88 | 0.214 | 0.123--0.305 | 0.233 | 0.119--0.348 | 57.4% |

## Interpretation boundary

The publisher field combines the simulation-derived basin adjustment and spatial nonergodic site term. The empirical field is a station residual after the MF2013 AVS30 and D1400 terms and still contains average-path structure. Their association tests independent spatial concordance; it does not equate the two amplitudes or validate a production nonergodic GMM. The simulation uses crustal scenarios and one 3D velocity model, so the source-matched result is reported first and its spatial-block uncertainty is retained.
