# Nearby-station physical audit for the 10,467-mesh public-parameter run

- Meshes linked: 10,467
- Candidate nearby stations with SA(3.0 s) site residuals and site fields: 2,261
- Median nearest-station distance: 7.8 km
- Meshes with nearest station within 50 km: 92.8%

## Key contrasts

- Northeast-Hokkaido meshes have a median nearest-station D1400 of 212 m, compared with 75 m outside that band.
- Northeast-Hokkaido meshes have a median nearest-station SA(3.0 s) site residual of -0.144, compared with -0.201 outside that band.
- Residual-tail meshes have a median nearest-station D1400 of 165 m, compared with 122 m for non-tail meshes.
- Residual-tail meshes have a median nearest-station SA(3.0 s) site residual of -0.191, compared with -0.194 for non-tail meshes.
- High-official-PBV meshes have a tail fraction of 26.6%.

## Interpretation

The residual-tail structure is geographically coherent and is colocated with nearby station fields that carry long-period site information. This supports a physical-attribution check based on deep-basin and station-residual proxies, while keeping the claim at the level of an audit. The calculation does not infer full 3-D basin geometry at every 250 m mesh.

## Output files

- `outputs/jshis_public_psha_10467_near_station_physical_link.csv`
- `outputs/jshis_public_psha_10467_near_station_physical_summary.csv`
