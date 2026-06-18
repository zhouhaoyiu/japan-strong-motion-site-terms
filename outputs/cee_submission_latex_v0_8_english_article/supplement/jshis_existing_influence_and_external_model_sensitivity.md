# Existing influence and external-model sensitivity summary

Generated on 2026-06-14.

## Event-exclusion sensitivity

| scenario_id | scenario_label | target | target_label | excluded_event_count | excluded_event_ids | basic_records | site_records | basic_mae | site_mae | mae_reduction | percent_mae_reduction | baseline_percent_mae_reduction | delta_vs_baseline_pct_points |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | All events | maxaccrd050 | PGA RotD50 | 0 |  | 322061 | 321732 | 0.4088 | 0.4022 | 0.006598 | 1.614 | 1.614 | 0 |
| all | All events | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 0 |  | 322061 | 321732 | 0.4156 | 0.3995 | 0.01612 | 3.879 | 3.879 | 0 |
| all | All events | rsaccrd050d005t0100 | SA(1.0s) RotD50 | 0 |  | 322061 | 321732 | 0.3855 | 0.3235 | 0.06194 | 16.07 | 16.07 | 0 |
| all | All events | rsaccrd050d005t0300 | SA(3.0s) RotD50 | 0 |  | 322061 | 321732 | 0.3928 | 0.3096 | 0.08315 | 21.17 | 21.17 | 0 |
| exclude_2011_march_sequence | Exclude 2011-03-11 to 2011-03-31 events | maxaccrd050 | PGA RotD50 | 46 | 35499;35504;35528;35536;35545;35848;36138;36177;36195;36254;36549;36647;36686;36727;36871;36897;37194;37239;37249;37507;37657;37665;37744;37752;37839;37898;37902;38015;38107;38168;... | 305418 | 305089 | 0.4086 | 0.4021 | 0.006546 | 1.602 | 1.614 | -0.01179 |
| exclude_2011_march_sequence | Exclude 2011-03-11 to 2011-03-31 events | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 46 | 35499;35504;35528;35536;35545;35848;36138;36177;36195;36254;36549;36647;36686;36727;36871;36897;37194;37239;37249;37507;37657;37665;37744;37752;37839;37898;37902;38015;38107;38168;... | 305418 | 305089 | 0.4146 | 0.3986 | 0.01606 | 3.873 | 3.879 | -0.006695 |

## Station-exclusion sensitivity

| scenario_id | scenario_label | target | target_label | excluded_site_count | excluded_siteids | basic_records | site_records | basic_mae | site_mae | mae_reduction | percent_mae_reduction | baseline_percent_mae_reduction | delta_vs_baseline_pct_points |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | All sites | maxaccrd050 | PGA RotD50 | 0 |  | 322061 | 321732 | 0.4088 | 0.4022 | 0.006598 | 1.614 | 1.614 | 0 |
| all | All sites | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 0 |  | 322061 | 321732 | 0.4156 | 0.3995 | 0.01612 | 3.879 | 3.879 | 0 |
| all | All sites | rsaccrd050d005t0100 | SA(1.0s) RotD50 | 0 |  | 322061 | 321732 | 0.3855 | 0.3235 | 0.06194 | 16.07 | 16.07 | 0 |
| all | All sites | rsaccrd050d005t0300 | SA(3.0s) RotD50 | 0 |  | 322061 | 321732 | 0.3928 | 0.3096 | 0.08315 | 21.17 | 21.17 | 0 |
| exclude_top1 | Exclude top 1 sites by record count | maxaccrd050 | PGA RotD50 | 1 | 1208041 | 321414 | 321085 | 0.4092 | 0.4026 | 0.006633 | 1.621 | 1.614 | 0.007108 |
| exclude_top1 | Exclude top 1 sites by record count | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 1 | 1208041 | 321414 | 321085 | 0.416 | 0.3998 | 0.01618 | 3.888 | 3.879 | 0.008833 |

## Magnitude-bin sensitivity

| mw_bin_id | mw_bin_label | mw_min_inclusive | mw_max_exclusive | target | target_label | basic_events | site_events | basic_records | site_records | basic_mae | site_mae | mae_reduction | percent_mae_reduction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mw_lt_5_0 | Mw<5.0 | 0 | 5 | maxaccrd050 | PGA RotD50 | 261 | 261 | 44414 | 44364 | 0.4239 | 0.4183 | 0.005581 | 1.317 |
| mw_lt_5_0 | Mw<5.0 | 0 | 5 | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 261 | 261 | 44414 | 44364 | 0.4175 | 0.4063 | 0.01118 | 2.678 |
| mw_lt_5_0 | Mw<5.0 | 0 | 5 | rsaccrd050d005t0100 | SA(1.0s) RotD50 | 261 | 261 | 44414 | 44364 | 0.3842 | 0.3179 | 0.06628 | 17.25 |
| mw_lt_5_0 | Mw<5.0 | 0 | 5 | rsaccrd050d005t0300 | SA(3.0s) RotD50 | 261 | 261 | 44414 | 44364 | 0.4256 | 0.3366 | 0.089 | 20.91 |
| mw_5_0_5_5 | 5.0<=Mw<5.5 | 5 | 5.5 | maxaccrd050 | PGA RotD50 | 845 | 845 | 137686 | 137570 | 0.416 | 0.4094 | 0.006547 | 1.574 |
| mw_5_0_5_5 | 5.0<=Mw<5.5 | 5 | 5.5 | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 845 | 845 | 137686 | 137570 | 0.4127 | 0.3971 | 0.01562 | 3.786 |

## Distance-bin sensitivity

| distance_bin_id | distance_bin_label | distance_min_km | distance_max_km | target | target_label | basic_records | site_records | basic_mae | site_mae | mae_reduction | percent_mae_reduction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| r_0_50 | 0<fault_dist<=50 km | 0 | 50 | maxaccrd050 | PGA RotD50 | 6513 | 6479 | 0.3982 | 0.4005 | -0.002287 | -0.5743 |
| r_0_50 | 0<fault_dist<=50 km | 0 | 50 | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 6513 | 6479 | 0.4542 | 0.4421 | 0.01209 | 2.662 |
| r_0_50 | 0<fault_dist<=50 km | 0 | 50 | rsaccrd050d005t0100 | SA(1.0s) RotD50 | 6513 | 6479 | 0.4541 | 0.3793 | 0.07484 | 16.48 |
| r_0_50 | 0<fault_dist<=50 km | 0 | 50 | rsaccrd050d005t0300 | SA(3.0s) RotD50 | 6513 | 6479 | 0.3704 | 0.2841 | 0.08635 | 23.31 |
| r_50_100 | 50<fault_dist<=100 km | 50 | 100 | maxaccrd050 | PGA RotD50 | 37633 | 37566 | 0.4151 | 0.4164 | -0.001279 | -0.3082 |
| r_50_100 | 50<fault_dist<=100 km | 50 | 100 | rsaccrd050d005t0030 | SA(0.3s) RotD50 | 37633 | 37566 | 0.4324 | 0.4218 | 0.01055 | 2.44 |

## Zhao et al. (2006) external GMPE overall metrics

| target | target_label | period | model_family | n_records | residual_sum | residual_sq_sum | absolute_error_sum | observed_sum | predicted_sum | mean_residual_ln | mae_ln | rmse_ln | mean_observed_ln_g | mean_predicted_ln_g | residual_std_ln | mae_log10_equiv | rmse_log10_equiv |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| maxaccrd050 | PGA RotD50 | PGA | zhao2006_openquake | 322020 | -1.891e+05 | 5.201e+05 | 3.264e+05 | -1.884e+06 | -1.695e+06 | -0.5871 | 1.014 | 1.271 | -5.852 | -5.265 | 1.127 | 0.4402 | 0.5519 |
| rsaccrd050d005t0030 | SA(0.3s) RotD50 | 0.30 | zhao2006_openquake | 322020 | -2.254e+05 | 5.545e+05 | 3.367e+05 | -1.678e+06 | -1.452e+06 | -0.6998 | 1.045 | 1.312 | -5.21 | -4.51 | 1.11 | 0.454 | 0.5699 |
| rsaccrd050d005t0100 | SA(1.0s) RotD50 | 1.00 | zhao2006_openquake | 322020 | -1.349e+05 | 3.017e+05 | 2.503e+05 | -2.045e+06 | -1.91e+06 | -0.419 | 0.7772 | 0.9679 | -6.349 | -5.93 | 0.8726 | 0.3375 | 0.4204 |
| rsaccrd050d005t0300 | SA(3.0s) RotD50 | 3.00 | zhao2006_openquake | 322020 | -1.423e+05 | 2.528e+05 | 2.324e+05 | -2.641e+06 | -2.499e+06 | -0.442 | 0.7216 | 0.886 | -8.202 | -7.76 | 0.7679 | 0.3134 | 0.3848 |

## Zhao et al. (2006) station-site correlations

| target | target_label | model_family | feature | n_sites | min_records_per_site | spearman_rho | pvalue |
| --- | --- | --- | --- | --- | --- | --- | --- |
| maxaccrd050 | PGA RotD50 | zhao2006_openquake | log10(VS10) | 2163 | 20 | -0.2951 | 9.886e-45 |
| maxaccrd050 | PGA RotD50 | zhao2006_openquake | log10(VS20) | 1669 | 20 | -0.3389 | 4.015e-46 |
| maxaccrd050 | PGA RotD50 | zhao2006_openquake | log10(VS30) | 1152 | 20 | -0.06118 | 0.03788 |
| maxaccrd050 | PGA RotD50 | zhao2006_openquake | log10(AVS30) | 2267 | 20 | -0.1919 | 3.094e-20 |
| maxaccrd050 | PGA RotD50 | zhao2006_openquake | log10(VS30_proxy) | 2267 | 20 | -0.2508 | 7.309e-34 |
| maxaccrd050 | PGA RotD50 | zhao2006_openquake | log10(D1100) | 2218 | 20 | -0.01535 | 0.4698 |

The current environment reports an OpenQuake geospatial-library import error, so no new external GMPE is generated here. The existing Zhao et al. (2006) derived outputs remain an external-model audit with explicit distance and VS30-proxy limitations.
