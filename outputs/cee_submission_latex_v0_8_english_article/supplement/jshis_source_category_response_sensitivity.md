# J-SHIS source-category response-spectrum sensitivity

## Data boundary

- The three inputs are official J-SHIS 2020 engineering-bedrock response-spectrum maps for all earthquakes, active shallow earthquakes and subduction earthquakes.
- all_earthquakes: `P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip`; SHA-256 `91ea26bdde63275f45c15e8a2a84d5382ca67d537468bedf494384700534862f`; https://www.j-shis.bosai.go.jp/map/respmap/data/P/Y2020/P-Y2020-RESP-MAP-AVR-TTL_MTTL-T50.zip
- active_shallow: `P-Y2020-RESP-MAP-AVR-LND_MTTL-T50.zip`; SHA-256 `ec1ab8436bffa628788a67b4cfeeded59b6debd9eb8979cdc82f95a0deeb1532`; https://www.j-shis.bosai.go.jp/map/respmap/data/P/Y2020/P-Y2020-RESP-MAP-AVR-LND_MTTL-T50.zip
- subduction: `P-Y2020-RESP-MAP-AVR-PPE_MTTL-T50.zip`; SHA-256 `75e9c127556a4b76da655c733f33cce0fbbb93a451f41d2b2513ffb6769d64f7`; https://www.j-shis.bosai.go.jp/map/respmap/data/P/Y2020/P-Y2020-RESP-MAP-AVR-PPE_MTTL-T50.zip
- Source-category ordinates at a fixed exceedance probability are evaluated separately. They are not summed, because probability aggregation is nonlinear.
- The same frozen spatial-block station multiplier is applied to each product. This isolates sensitivity to the official source category; it does not estimate source-specific station terms.

## SA(3.0 s), 10% in 50 years

- Paired stations: 1,628.
- Active-shallow ordinate is larger at 200 stations; subduction ordinate is at least as large at 1,428 stations.
- active_shallow: median official ordinate 0.031 g; median surface reference 0.032 g; median station-adjusted value 0.029 g.
- subduction: median official ordinate 0.071 g; median surface reference 0.072 g; median station-adjusted value 0.067 g.

## Interpretation

The calculation demonstrates the station correction on official maps spanning distinct source regimes. It remains a matched-station response-spectrum sensitivity. A production J-SHIS PSHA rerun would require the complete official source occurrence model, logic-tree weights and implementation details inside the hazard integral.
