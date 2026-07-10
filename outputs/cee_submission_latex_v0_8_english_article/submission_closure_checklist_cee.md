# CEE submission checklist

Date: 11 July 2026

## Manuscript checks

- [x] Title has 13 words and no punctuation.
- [x] Abstract has 140 words, contains no citations and follows the journal's background-to-"Here we present" structure.
- [x] Section order is Introduction, Results, Discussion and Methods.
- [x] Introduction and Discussion contain no subheadings.
- [x] Main text is below the 5,000-word limit.
- [x] Six main figures and one main table are cited in order.
- [x] Supplementary items are cited in numerical order.
- [x] Full postal address, correspondence email and ORCIDs are present.
- [x] Funding, competing interests, author contributions and acknowledgements are present.
- [x] Generative-AI assistance is disclosed in Methods and no generated image is included.
- [x] All 33 main references and 11 Supplementary Information references have a DOI or official HTTPS source and appear in first-citation order.
- [x] English and Chinese reference lists use the same 33 records in the same order.

## Scientific checks

- [x] Eight periods are calculated directly from flatfile columns.
- [x] Primary residual and spectrum analyses use official ground-surface installation records only.
- [x] MF2013 uses model-compatible RotD100; Zhao 2006 uses RotD50 as the closest available orientation-independent coordinate to its geometric-mean target.
- [x] RotD50--RotD100 station-field sensitivity is quantified on the same 222,664-record sample.
- [x] Event and station terms are separated and constrained to weighted zero means.
- [x] All decompositions converge below `1e-10` log10 units.
- [x] Event-group repeatability uses separate converged two-way fixed-effect fits on disjoint earthquake groups and a common event partition across periods and backbones.
- [x] Zhao 2006 repeats the decomposition, event holdout and spatial validation on the same surface observations.
- [x] Cross-network, fixed-macroregion and top-event deletion tests are reported with adverse results retained.
- [x] Path and source strata are solved by connected component, with solver residuals and event split-half repeatability reported.
- [x] Equal-stratum prediction retains negative spatial-fold results and is labelled as a sensitivity analysis.
- [x] Empirical station-model prediction intervals are calibrated outside each held spatial block.
- [x] Hazard propagation uses out-of-fold zero-centred predictions.
- [x] Official `Vs=400 m/s` values are converted to station AVS30 conditions before station adjustment.
- [x] No claim of a national directional spectrum shift remains.
- [x] No claim of an official source-level J-SHIS PSHA reproduction remains.

## Files

- [x] `main.pdf` compiled and visually checked.
- [x] `supplementary_information.pdf` compiled and visually checked.
- [x] Supplementary Information contains seven figures and twelve tables in first-citation order.
- [x] Derived tables pass the updated `work/validate_event_adjusted_release.py`.
- [x] Main and supplementary figures have source scripts.
- [x] Historical figures and obsolete PBV audit tables are removed from the release package.
- [x] `reference_traceability.md` records DOI or official-source checks and local full-text coverage.

## Author confirmations before upload

- [ ] Confirm final author names, order, affiliation, email and ORCIDs.
- [ ] Confirm all authors approve submission to *Communications Earth & Environment*.
- [ ] Confirm the GitHub repository is public and contains the same release.
- [ ] Mint and insert the Zenodo DOI if the authors choose to archive before first submission.
- [ ] Confirm the final data, code and AI-use statements in the submission system.
- [ ] Select article type, subject terms, reviewers and exclusions in the portal.
