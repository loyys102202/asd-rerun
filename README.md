# ASD blood expression score analysis

This repository contains the analysis code and workflow for the revised iScience manuscript:

> **A monocyte-tracking blood expression score reveals platform-specific dissociation between cellular composition and autism diagnosis** &mdash; Jing Wen (manuscript ID ISCIENCE-D-26-03472)

The study reanalyses public GEO pediatric blood transcriptome datasets to evaluate the transportability of an exploratory blood expression score across platforms, blood-fraction preparations, and analytic cohorts. The score is not intended as a clinically deployable diagnostic test.

## Scope and framing

This is a reproducibility and methods-transparency package supporting an **exploratory biomarker/signature discovery** study. The reported score:

- was trained in a single male-only pediatric whole-blood cohort (GSE18123-GPL570, n = 99),
- shows partial transportability to a same-series cross-platform cohort (GSE18123-GPL6244),
- does **not** discriminate autism in a fully independent leukocyte cohort (GSE42133),
- does **not** transport without re-fitting to a lymphocyte-enriched cohort (GSE25507),
- but tracks monocyte-marker expression consistently in all four cohorts.

The repository should not be interpreted as supporting a validated clinical classifier.

## Contents

- `01_download_geo.py` &mdash; download and unify the four GEO series
- `02_preprocess.py` &mdash; probe-to-gene IQR-max aggregation, log2 retention
- `03_common_genes.py` &mdash; cross-platform gene intersection (GPL570 &cap; GPL6244)
- `04_nested_cv.py` &mdash; repeated nested cross-validation (100 repeats &times; 5 outer folds)
- `05_lock_panel.py` &mdash; stability-selected locked 7-gene panel
- `06_bootstrap.py` &mdash; 0.632 bootstrap optimism correction (B = 1,000)
- `07_permutation_train.py` &mdash; label-permutation null on training cohort
- `08_external_GPL6244.py` &mdash; same-series cross-platform external test
- `09_external_GSE42133.py` &mdash; independent external test (training-cohort-refitted 6-gene reduced score, because ERVK3-2 is unavailable on Illumina HT-12 v4)
- `10_tissue_boundary_GSE25507.py` &mdash; tissue-boundary analysis on lymphocyte-enriched samples
- `11_cell_markers.py` &mdash; blood-cell marker correlation analysis (CIBERSORT LM22 + xCell reference panels)
- `12_*` and `figures/`, `tables/` &mdash; figure and table generation scripts

## Data

All datasets analysed in this study are public NCBI GEO datasets:

- GSE18123 (sub-cohorts on GPL570 and GPL6244)
- GSE42133
- GSE25507

No new patient samples were collected. No human-subject data are deposited in this repository.

## Code availability

The complete analysis workflow is provided in this repository. A permanent archived version is available at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20404078.svg)](https://doi.org/10.5281/zenodo.20404078)

Citation file: `CITATION.cff`. Random seed `42` is used throughout. End-to-end reproduction: `bash run_all.sh`.

## Environment

Python 3.11 with the package versions pinned in `requirements.txt` (scikit-learn 1.5, statsmodels 0.14, NumPy 1.26, pandas 2.2, SciPy 1.13, matplotlib 3.9).

## Important note

This repository supports an exploratory blood-transcriptome signature-discovery study. The reported score should not be interpreted as a validated clinical diagnostic classifier. Female generalisability is not established because the discovery cohort was male-only and female external testing was underpowered.

## License

MIT, see `LICENSE`.
