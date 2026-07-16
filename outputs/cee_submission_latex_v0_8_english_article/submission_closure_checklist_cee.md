# CEE submission checklist

Date: 11 July 2026

## Manuscript checks

- [x] Title has 13 words and no punctuation.
- [x] Abstract is below 150 words, contains no citations and moves directly from two background sentences to the analysis.
- [x] Section order is Introduction, Results, Discussion and Methods.
- [x] Introduction and Discussion contain no subheadings.
- [x] Main text is below the 5,000-word limit.
- [x] Eight main figures and one main table are cited in order.
- [x] All Supplementary Figures and Tables are cited and their numbering matches the compiled Supplementary Information.
- [x] Full postal address, correspondence email and ORCIDs are present.
- [x] Funding, competing interests, author contributions and acknowledgements are present.
- [x] Author contributions use full names, and references with more than five authors use first author plus `et al.`.
- [x] Generative-AI assistance is disclosed in Methods and no generated image is included.
- [x] All 37 main references and 16 Supplementary Information references have a DOI or official HTTPS source and appear in first-citation order.
- [x] English and Chinese reference lists use the same 37 records in the same order.
- [x] Current CEE guidance was rechecked on 11 July 2026: the initial submission uses compiled PDFs with figures embedded, main text is below the recommended 5,000 words excluding Methods, and LLM use is documented in Methods.

## Scientific checks

- [x] Eight periods are calculated directly from flatfile columns.
- [x] The `sub1-v2024` JMA-magnitude threshold and the later MF2013 finite-`Mw` screen are distinguished, with the full 333,808-to-222,664 record flow audited.
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
- [x] Crustal, interplate and intraplate station fields use a fixed 20-record threshold, unchanged spatial models and all eight periods; adverse crustal and interplate periods are retained.
- [x] The Sung et al. Kanto comparison uses the publisher-supplied field, fixed nearest-grid matching, natural-log unit conversion and 31.7 km spatial-block intervals; the weak source-conditioned out-of-fold result is reported.
- [x] Empirical station-model prediction intervals are calibrated outside each held spatial block.
- [x] Six fixed spatial-model complexity settings retain positive overall gains at all eight periods and are reported as sensitivity tests without model reselection.
- [x] The full Kanno 2006 shallow/deep PGA cutoff is reconstructed and included with the MF2013 magnitude, distance and event-station screens; the adverse 0.1 s result and sample-dependent aggregate median are stated explicitly.
- [x] The national KiK-net paired-sensor validation uses 102,428 surface--borehole spectra and disjoint-event calibration and evaluation.
- [x] Cross-period paired-sensor spectral shape is validated with station-cluster bootstrap intervals and disjoint event groups.
- [x] Frozen K-NET models are evaluated on KiK-net station terms from disjoint earthquakes, with no target recentering or target labels used in training.
- [x] Chronological 70/30, 80/20 and 90/10 event splits are reported in both transfer directions, with an event-level split manifest and the adverse reverse result retained.
- [x] Post-2015 ESM records reproduce 1--3 s station recurrence within the published Bindi magnitude range, with station-bootstrap intervals and an independent Kotha/RotD50 sensitivity.
- [x] The adverse Japanese-to-ESM transfer using harmonised VS30 and elevation is retained and interpreted as a regional prediction boundary.
- [x] Local ESM spatial-block and held-event-plus-spatial-block tests use a response-independent minimum-block rule; all adverse 1--3 s aggregate results are retained.
- [x] Independent correction-field agreement and empirical prediction intervals are reported as distinct uncertainty quantities.
- [x] Hazard propagation uses out-of-fold zero-centred predictions.
- [x] Official `Vs=400 m/s` values are converted to station AVS30 conditions before station adjustment.
- [x] Official active-shallow and subduction response-spectrum products are evaluated separately; fixed-probability category ordinates are not summed.
- [x] Source-conditioned crustal, interplate and intraplate fields are propagated separately and are not combined without official source weights.
- [x] No claim of a national directional spectrum shift remains.
- [x] No claim of an official source-level J-SHIS PSHA reproduction remains.

## Files

- [x] `main.pdf` compiled and visually checked.
- [x] `supplementary_information.pdf` compiled and visually checked.
- [x] Supplementary Information contains fourteen figures and twenty-six tables.
- [x] Derived tables pass the updated `work/validate_event_adjusted_release.py`.
- [x] Main and supplementary figures have source scripts.
- [x] Historical figures and obsolete PBV audit tables are removed from the release package.
- [x] `reference_traceability.md` records DOI or official-source checks and local full-text coverage.
- [x] A compact peer-review archive provides code, the input manifest and derived summary tables with SHA-256 checksums.

## Author confirmations before upload

- [ ] Confirm final author names, order, affiliation, email and ORCIDs.
- [ ] Confirm all authors approve submission to *Communications Earth & Environment*.
- [ ] Make the reviewed GitHub release public and archive it with a DOI before publication.
- [ ] Confirm the final data, code and AI-use statements in the submission system.
- [ ] Select article type, subject terms, reviewers and exclusions in the portal.
