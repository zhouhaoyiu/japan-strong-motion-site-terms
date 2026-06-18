# Multi-event KiK-net transfer-function audit

## Data inventory

- Event-station pairs: 460.
- Events: 2.
- Stations: 430.
- Frequency-bin rows: 14,688.

## Event coverage

| event_id | source_formats | station-events | median horizontal PGA ratio | median 0.2-0.5 Hz ratio | median 1-2 Hz ratio |
| --- | --- | ---: | ---: | ---: | ---: |
| 2512082315 | converted_mseed | 288 | 1.111 | 2.382 | 1.152 |
| 2601061018 | nied_ascii | 172 | 2.352 | 1.053 | 1.522 |

## Combined transfer-function result

- Median horizontal PGA ratio: 1.458.
- Interquartile horizontal PGA ratio: 0.902 to 2.223.
- Event-station pairs with ratio > 1: 333 / 460.
- Event-station pairs with ratio > 2: 138 / 460.

## Boundary

This audit estimates observed surface/downhole transfer functions from two available KiK-net events. The 2026 ASCII event is converted with the NIED header scale factors. The 2025 converted miniSEED event is treated as a same-record amplitude-ratio data path, consistent with the existing local conversion. The result is a multi-event physical consistency check for station terms and frequency bands; it is not a nonlinear site-response model and does not replace a full KiK-net archive analysis.

## Output tables

- `outputs/kiknet_multievent_transfer_functions.csv`
- `outputs/kiknet_multievent_transfer_function_bins.csv`
- `outputs/kiknet_multievent_transfer_function_event_summary.csv`
- `outputs/kiknet_multievent_transfer_function_station_summary.csv`
