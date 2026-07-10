# J-SHIS flatfile selection audit

- Archive: `work/external_data/jshis_gmf/flatfile_sub1-v2024.zip`
- Record rows and unique `smrec_id` values: 333,808
- Missing site joins: 0
- Missing source joins: 0
- JMA magnitude range: 5.0--9.0
- Shortest-fault-distance range: 1.0001--299.9953 km
- Surface records with nonpositive AVS30: 42
- Nonpositive-AVS30 records that also lack finite Mw: 1

## Sequential selection

| Stage | Records | Excluded at stage | Earthquakes | Stations |
|---|---:|---:|---:|---:|
| Public sub1-v2024 archive | 333,808 | 0 | 1,840 | 2,581 |
| Ground-surface installation | 231,380 | 102,428 | 1,840 | 1,882 |
| Finite F-net moment magnitude Mw | 222,705 | 8,675 | 1,737 | 1,881 |
| Supported source class and positive distance | 222,705 | 0 | 1,737 | 1,881 |
| Positive AVS30 | 222,664 | 41 | 1,737 | 1,880 |
| Finite RotD100 and MF2013 predictions at all eight periods | 222,664 | 0 | 1,737 | 1,880 |

## Period-specific final counts

0.1 s: 222,664, 0.2 s: 222,664, 0.3 s: 222,664, 0.5 s: 222,664, 1 s: 222,664, 2 s: 222,664, 3 s: 222,664, 5 s: 222,664

Every record in the downloaded subset satisfies `mjma >= 5`. MF2013 uses the separate F-net `mw` field, which creates the magnitude-related loss between the surface archive and the model sample.
