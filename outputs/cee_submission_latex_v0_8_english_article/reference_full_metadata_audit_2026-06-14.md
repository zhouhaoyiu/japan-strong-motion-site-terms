# Full Reference Metadata Audit

Date: 2026-06-17.

Scope: references in `outputs/cee_submission_latex_v0_8_english_article/main.tex`.

## Result

The current CEE manuscript contains 14 numbered references. The list was rechecked against publisher pages, DOI metadata, J-SHIS/NIED official pages, OpenQuake documentation, GSI documentation, and the Nozu (2023) article page on 2026-06-17. Thirteen items have DOI metadata or official web-source metadata. The only remaining high-risk item is the WCEE 2024 flatfile proceedings reference, which is retained because it is the official J-SHIS/NIED flatfile-construction reference and no DOI or stable paper number was found during the latest check.

No duplicate de la Torre reference remains in the manuscript. The list now includes the added GSI elevation-tile, J-SHIS WMS surface-ground model, and Nozu (2023) references used by the terrain, site-amplification, and northern-Hokkaido mechanism audits.

## Item-by-item check

| No. | Working reference | Verification source | Verified fields | Action |
|---:|---|---|---|---|
| 1 | Morikawa and Fujiwara 2013 | Fuji Technology Press page, J-SHIS MF2013 page, DOI `10.20965/jdr.2013.p0878` | title, authors, journal, volume 8, pages 878--888, year 2013 | no change |
| 2 | Zhao et al. 2006 | BSSA/GeoScienceWorld and DOI `10.1785/0120050122` | title, journal, volume 96, issue 3, pages 898--913, year 2006 | no change; `et al.` is acceptable in numbered manuscript style |
| 3 | Stewart, Afshari and Goulet 2017 | SAGE/GeoScienceWorld page and DOI `10.1193/081716EQS135M` | title, authors, journal, volume 33, pages 1385--1414, year 2017 | no change |
| 4 | Lavrentiadis et al. 2023 | Springer page and DOI `10.1007/s10518-022-01485-x` | title, journal, volume 21, pages 5121--5150, year 2023 | no change |
| 5 | Lacour 2023 | Springer page and DOI `10.1007/s10518-022-01402-2` | title, journal, volume 21, pages 5209--5232, year 2023 | no change |
| 6 | NIED Strong Ground Motion Flat File 2025 | J-SHIS/NIED flatfile page and DOI `10.17598/NIED.0032` | official English title, publisher, year, DOI; page confirms 2024 version | no change |
| 7 | Morikawa et al. 2024 WCEE flatfile proceedings | J-SHIS/NIED flatfile page and web search | author lead, title, WCEE container, Milan, 30 June--5 July 2024; no DOI/paper number found | retain; final reference-manager export recommended |
| 8 | J-SHIS/NIED MF2013 page | J-SHIS/NIED MF2013 official page | page title, URL, year/access basis | no change |
| 9 | Pagani et al. 2014 | SRL/GeoScienceWorld, OpenQuake bibliography, DOI `10.1785/0220130087` | title, authors, journal, volume 85, issue 3, pages 692--702, year 2014 | no change |
| 10 | Parker and Stewart 2022 | Earthquake Spectra publisher page and DOI `10.1177/87552930211056963` | title, authors, journal, volume 38, pages 841--864, year 2022 | no change |
| 11 | de la Torre et al. 2024 | SAGE/Wiley page and DOI `10.1177/87552930241270562` | title, author list, journal, volume 40, pages 2475--2503, year 2024 | no change |
| 12 | GSI elevation tile specification | GSI Maps elevation-tile specification page | official page title and URL for elevation-tile source | no change; web source has no DOI |
| 13 | J-SHIS WMS public surface-ground model fields | J-SHIS WMS URL used by the site-amplification audit | official model-field URL and access date | no change; web source has no DOI |
| 14 | Nozu 2023 | Earth, Planets and Space article page and DOI `10.1186/s40623-023-01936-y` | title, author, journal, volume 75, article 177, year 2023, DOI | no change |

## Notes for final submission

- The English CEE manuscript uses numbered Nature-style references. DOI values are printed only for data/web items where the current manuscript style already uses them.
- Ref. 7 remains the only reference requiring final bibliographic confirmation. If the WCEE paper ID becomes available, replace the current proceedings entry with the official conference export.
- Web-source access dates should be refreshed on the actual submission date if the target journal requests access dates.

## Sources checked

- Fuji Technology Press JDR page for Morikawa and Fujiwara (2013): `https://www.fujipress.jp/jdr/dr/dsstr000800050878/`
- J-SHIS/NIED Strong Ground Motion Flat File page: `https://www.j-shis.bosai.go.jp/en/labs/ground-motion-flatfile/`
- J-SHIS/NIED MF2013 page: `https://www.j-shis.bosai.go.jp/en/labs/mf2013/`
- Stewart et al. (2017) publisher page: `https://journals.sagepub.com/doi/10.1193/081716eqs135m`
- Lavrentiadis et al. (2023) publisher page: `https://link.springer.com/article/10.1007/s10518-022-01485-x`
- Lacour (2023) publisher page: `https://link.springer.com/article/10.1007/s10518-022-01402-2`
- OpenQuake bibliography for Pagani et al. (2014): `https://docs.openquake.org/oq-engine/3.21/manual/underlying-science/bibliography.html`
- Parker and Stewart (2022) publisher page: `https://pubs.geoscienceworld.org/eeri/earthquake-spectra/article/38/2/841/613712/Ergodic-site-response-model-for-subduction-zone`
- de la Torre et al. (2024) publisher page: `https://journals.sagepub.com/doi/10.1177/87552930241270562`
- GSI elevation-tile specification: `https://maps.gsi.go.jp/development/demtile.html`
- Nozu (2023) article page: `https://earth-planets-space.springeropen.com/articles/10.1186/s40623-023-01936-y`
