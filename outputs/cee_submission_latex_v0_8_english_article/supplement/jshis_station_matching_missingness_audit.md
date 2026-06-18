# Station matching and missingness audit

Generated on 2026-06-14. Counts come from derived public-data audit tables.

| check | count | denominator | percent | source | note |
| --- | --- | --- | --- | --- | --- |
| strong_motion_records | 333808 |  |  | jshis_smrec_schema_sub1_summary.md | Public flatfile records in sub1-v2024. |
| unique_record_sites | 2581 |  |  | jshis_smrec_schema_sub1_summary.md | Unique site_id//10 values in the strong-motion records. |
| site_schema_rows | 2607 |  |  | jshis_site_schema_v2024_sub1.csv | Rows in the public site schema. |
| record_site_schema_matched_rows | 333808 | 3.338e+05 | 100 | jshis_smrec_schema_sub1_summary.md | Record-to-site match reported by the schema extraction audit. |

The CSV file contains variable-level missingness, MF2013 residual coverage, minimum-record thresholds, and official SA(3.0 s) response-map matching.
