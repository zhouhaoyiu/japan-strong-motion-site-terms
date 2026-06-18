# Official J-SHIS PSHA Input-Chain Audit

Audit date: 2026-06-14. Official model version: Y2024.

## What Was Checked

The audit downloads the official J-SHIS Probabilistic Seismic Hazard Maps parameter packages and records their file roles, field structures, and reconstruction status. The raw J-SHIS zip files stay in the local `work/` cache and are not redistributed in the manuscript package.

- CSV parameter package: `P-Y2024-PRM.zip` from https://www.j-shis.bosai.go.jp/map/JSHIS2/data/P/Y2024/P-Y2024-PRM.zip
- Shapefile parameter package: `P-Y2024-PRM-SHAPE.zip` from https://www.j-shis.bosai.go.jp/map/JSHIS2/data/P/Y2024/PRM/P-Y2024-PRM-SHAPE.zip
- J-SHIS download page: https://www.j-shis.bosai.go.jp/map/JSHIS2/download.html?lang=en
- J-SHIS API list: https://www.j-shis.bosai.go.jp/en/api-list
- File-format specification: https://www.j-shis.bosai.go.jp/map/JSHIS2/data/DOC/DataFileRule/A-RULES_en.pdf

## Main Audit Result

The official 2024 package contains 421 CSV parameter files (78.75 MB) and 200 shapefile layers (800 component files, 20.06 MB). The CSV package includes source activity parameters, source geometries, occurrence-frequency grids for sources without specified faults, zoning polygons, interplate/intraplate ratios, attenuation-relation parameter codes, and Pacific/Philippine Sea plate geometry.

This is stronger evidence than using only J-SHIS hazard-curve API outputs. It confirms that major official PSHA input layers are publicly accessible and versioned.

The audit still does not constitute an official full J-SHIS PSHA rerun. A full rerun would require converting these J-SHIS-specific parameter files into an executable PSHA source model, validating the attenuation-relation implementation and calculation switches, and reproducing the official probability-case aggregation for all national grid cells. The manuscript should describe the present addition as an official input-chain audit and bounded reconstruction check, not as a completed official national hazard-map recomputation.

## Category Summary

| category | n_csv_files | example_file | fields_or_block_structure | reconstruction_status | paper_use |
| --- | --- | --- | --- | --- | --- |
| source_geometry | 237 | P-Y2024-PRM-SHP_TYPE1_LND_A98F.csv | Source code, magnitude or frequency block, point/fault-plane coordinates, depth, length, width, strike, dip | Public input in CSV and shapefile forms. Geometry still needs source-model conversion for a full rerun. | Shows that official source geometry is downloadable and versioned. |
| activity_background | 107 | P-Y2024-PRM-ACT_LND_CGR5_CRUST_CV_SM.csv | MNO, JLG, JLA, WLG, WLA, FRQ, BVL, MMN, ANO, DEP, STR, DIP | Public input. Needs gridded-source translation and validation against J-SHIS calculation settings. | Documents the public background-source input layer. |
| zoning_area | 65 | P-Y2024-PRM-AREA_SHP_LND_CGR5_01.csv | JLON, JLAT, WLON, WLAT | Public input. Needed to reconstruct area-source occurrence models. | Documents spatial zoning used by background-source files. |
| activity_characteristic | 5 | P-Y2024-PRM-ACT_AVR_LND_A98F.csv | CODE, PROC, AVRACT, NEWACT, ALPHA, P_T30, P_T50, NAME | Public input. Needs conversion from J-SHIS CSV conventions to an executable PSHA source model. | Supports a documented official-source input-chain audit. |
| eqthr | 2 | P-Y2024-PRM-AVR_LND_A98F_EQTHR.csv | Fault code, recurrence interval, magnitude bounds, b-value, and fault geometry | Public input. Needs conversion and validation. | Documents the active-fault supplement layer. |
| inter_intra_ratio | 2 | P-Y2024-PRM-RATIO_INTER_INTRA.csv | EQCODE, ANO, INTERR, INTRAR | Public input. Must be mapped into the executable source logic. | Documents plate-source component ratios. |
| plate_shape | 2 | P-Y2024-PRM-PLATE_SHP_PSE_CPCF.csv | MNO, JLG, JLA, WLG, WLA, DEP | Public input. Needs validated interpolation and source-placement logic. | Documents public plate-geometry inputs. |
| attenuation_formula | 1 | P-Y2024-PRM-ATTENUATION_FORMULA.csv | EQCODE, EQTYPE, SPTYPE, MTTYPE, CRTYPE | Public input. The numerical GMPE implementation and calculation switches must be verified before claiming a full official rerun. | Defines the remaining bridge between official source inputs and ground-motion calculation. |

## Download Manifest

| asset | url | local_cache_path | redistributed | use_in_this_study |
| --- | --- | --- | --- | --- |
| P-Y2024-PRM.zip | https://www.j-shis.bosai.go.jp/map/JSHIS2/data/P/Y2024/P-Y2024-PRM.zip | work/jshis_full_psha_audit/downloaded/P-Y2024-PRM.zip | no | Local official parameter audit only; derived inventory tables are redistributed. |
| P-Y2024-PRM-SHAPE.zip | https://www.j-shis.bosai.go.jp/map/JSHIS2/data/P/Y2024/PRM/P-Y2024-PRM-SHAPE.zip | work/jshis_full_psha_audit/downloaded/P-Y2024-PRM-SHAPE.zip | no | Local official geometry audit only; derived inventory tables are redistributed. |
| J-SHIS file format specification | https://www.j-shis.bosai.go.jp/map/JSHIS2/data/DOC/DataFileRule/A-RULES_en.pdf | work/jshis_full_psha_audit/docs/A-RULES_en.pdf | no | Field definitions and file naming rules. |

## Recommended Manuscript Wording

We audited the official J-SHIS 2024 PSHM parameter packages to check whether the station-residual sensitivity analysis could be connected to the official PSHA input chain. The packages provide downloadable source activity parameters, source-geometry files, occurrence-frequency grids for earthquakes without specified source faults, zoning polygons, interplate/intraplate ratios, attenuation-relation parameter codes, and Pacific/Philippine Sea plate geometry. This audit supports the provenance of the official hazard-facing check and identifies the remaining step needed for a full official rerun: translation of the J-SHIS parameter conventions into a validated executable PSHA model with the official calculation settings. The present manuscript keeps the hazard calculation at the level directly supported by public products: official response-spectrum ordinates, official hazard-curve API references, official parameter-chain provenance, and an independent OpenQuake/GEM Japan-model sensitivity run.

## Files Written

- `jshis_official_psha_input_audit_summary.csv`
- `jshis_official_psha_input_file_inventory.csv`
- `jshis_official_psha_input_shape_inventory.csv`
- `jshis_official_psha_input_download_manifest.csv`
- `figures/jshis_official_psha_input_chain_audit.png` and `.pdf`
