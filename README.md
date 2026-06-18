# Japan Strong-Motion Site-Term Manuscript Package

This repository is a curated pre-submission workspace for a research project on long-period site and basin residuals in Japanese strong-motion records and their sensitivity for uniform-hazard spectrum ordinates.

The upload intentionally excludes raw third-party data, local caches, virtual environments, staged rerun directories, large exploratory waveform tables, and private planning notes. The current public-facing materials are:

- `outputs/cee_submission_latex_v0_8_english_article/`: English CEE-style manuscript source, compiled PDFs, figures, supplement, cover letter draft, and derived tables.
- `outputs/cee_v08_public_archive_readiness.md`: pre-deposition public-archive readiness review for the current v0.8 CEE package.
- `outputs/cee_v08_public_archive_file_review.csv`: file-level release review table for scripts, figures, manuscript files, and derived tables.
- `outputs/cee_public_release_candidate_v0_8_audit.md`: Git-visible public release candidate audit; the report records the current file count and raw/cache, local absolute-path, private-planning text, and zip member-name scan results.
- `outputs/cee_public_release_candidate_v0_8_smoke_test.md`: smoke-test report for the local public release candidate; the report records the current gate count, including nested manuscript-zip checksum, PDF page counts, sensitive-text scan, and Python syntax compilation.
- `outputs/cee_public_release_candidate_v0_8_sha256.txt`: SHA-256 checksum for the local public release candidate zip, kept outside the zip to avoid a self-referential archive hash.
- `outputs/cee_release_metadata_zenodo_v0_8_draft.json`: Zenodo metadata draft to finalize after author, funding, license, and repository approval.
- `outputs/cee_submission_package_sha256.txt`: SHA-256 checksum for the current v0.8 manuscript archive.

## Data Scope

The J-SHIS/NIED strong-motion flatfile, full response-spectrum map zip files, and other raw third-party inputs are not committed to this repository. Use the documented source links and reproduction instructions in the manuscript and reproducibility package to obtain raw inputs from official sources.

The OpenQuake/GEM Japan-model calculation in this repository is an independent PSHA sensitivity check for representative stations. It should not be described as an official J-SHIS full hazard-map reproduction.

## Current Status

This is not a final journal submission package. The manuscript package now includes the official J-SHIS response-spectrum sensitivity check, spatial and bootstrap uncertainty checks, positive-correction station audit, additional validation experiments, an independent OpenQuake/GEM Japan-model PSHA sensitivity run, a 10,467-mesh public-parameter J-SHIS PBV execution audit, a near-station physical audit linking residual-tail meshes to public site proxies, a basin-period proxy audit for the 1--3 s residual-tail structure, a tail-robustness case-profile audit, a nearest-station official borehole-profile audit using public NIED station profiles, a sampled-grid SA(3.0 s) continuous correction-surface audit, an independent GSI terrain geomorphic proxy audit, a J-SHIS VS400 gridded site-amplification audit, and a v0.8 public-archive readiness review. Final author approval, funding/COI statements, and public DOI or repository citation remain required before formal submission.

## Research Integrity

Do not treat these materials as proof that the package is ready for upload. Claims should remain limited to the analyses documented in the package, with the stated assumptions and limitations preserved.
