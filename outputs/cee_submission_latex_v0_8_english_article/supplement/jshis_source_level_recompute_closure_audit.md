# Public J-SHIS source-parameter boundary check

This audit separates public source-parameter and sampled-grid consistency checks from a full official J-SHIS production rerun.

## Closure state

| item | value | status |
| --- | --- | --- |
| official PRM source files inventoried | 421 | closed |
| official shapefile layers inventoried | 200 | closed |
| national sampled meshes checked | 10467 | closed |
| official threshold rows checked | 41868 | closed |
| failed or zero-probability rows | 0 | closed |
| official spatial correlation and logic-tree implementation | not public in executable form | open |
| validated official GMPE switch parity | not public in executable form | open |

## Claim boundary

| claim_class | statement |
| --- | --- |
| Can be claimed | Public J-SHIS source-input inventory and sampled-grid public-parameter checks on 10,467 meshes. |
| Can be claimed | All-period official response-spectrum ordinate propagation using observed station terms. |
| Open boundary | Completed official J-SHIS national production PSHA rerun. |
| Reason | The public packages lack a validated executable implementation of every official calculation switch, spatial correlation, and logic-tree aggregation step. |

## Manuscript wording

The public source-parameter packages support an input-chain inventory and sampled-grid consistency check. The official response-spectrum sensitivity remains tied to official J-SHIS map ordinates, while source recurrence, spatial correlation, model uncertainty, and logic-tree aggregation remain in the official product.
