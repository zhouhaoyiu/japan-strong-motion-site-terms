# Spatial station-model complexity sensitivity

The six settings use the same station targets, spatial blocks, preprocessing, weights and random seed (20260710). Only tree complexity, regularization or iteration budget changes. The primary setting remains fixed; the alternatives are a sensitivity audit and are not used to select a replacement model.

| Setting | SA(3.0 s) RMSE gain (%) | SA(3.0 s) correlation |
|---|---:|---:|
| fewer_leaves | 13.98 | 0.512 |
| larger_minimum_leaf | 10.71 | 0.483 |
| more_leaves | 10.45 | 0.467 |
| primary | 12.56 | 0.500 |
| shorter_iteration_budget | 12.85 | 0.497 |
| stronger_l2 | 12.29 | 0.496 |

Across all periods and settings, overall RMSE gains range from 7.07% to 19.27%.
The primary setting reproduces eight period rows and its minimum gain is 10.07%.
