# asd-rerun v3 — full reproducibility code for ISCIENCE-D-26-03472

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20404078.svg)](https://doi.org/10.5281/zenodo.20404078)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

Complete pipeline + figures + tables for the manuscript
*"A monocyte-tracking blood expression score reveals platform-specific
dissociation between cellular composition and autism case-control status"*
(Jing Wen et al., *iScience*, in revision; manuscript ID **ISCIENCE-D-26-03472**).

> **Archived release**: this code is archived on Zenodo at
> [https://doi.org/10.5281/zenodo.20404078](https://doi.org/10.5281/zenodo.20404078)
> (version v1.0.0).

## Quick start

```bash
pip install -r requirements.txt

# Option A — full reproduction from GEO (downloads ~1 GB raw data):
bash run_all.sh

# Option B — just regenerate figures and tables from precomputed
# artifacts (extract asd-data.zip into this directory first):
bash run_all.sh figures
```

## Repository structure

```
asd-rerun/
├── README.md, LICENSE, requirements.txt, run_all.sh
│
├── config.py                          cohort definitions + hyperparameters
├── pipeline.py                        shared utilities (download, preprocess, scoring)
│
├── 01_download_geo.py                 fetch GSE18123 + GSE42133 + GSE25507
├── 02_preprocess.py                   log2 → IQR-max probe-to-gene
├── 03_common_genes.py                 GPL570 ∩ GPL6244 = 17,923 common genes
├── 04_nested_cv.py                    100 × 5-fold nested CV + stability selection
├── 05_lock_panel.py                   lock 7-gene panel at 60 % stability threshold
├── 06_bootstrap.py                    .632 bootstrap optimism correction
├── 07_permutation_train.py            training permutation null (1000 perms)
├── 08_external_GPL6244.py             platform-adapted + strict + recalibration
├── 09_external_GSE42133.py            6-gene refit + 7-gene zero-fill
├── 10_tissue_boundary_GSE25507.py     adapted + within-tissue refit
├── 11_cell_markers.py                 8 cell-type marker correlations
├── 11b_precompute_fold_change.py      bootstrap fold change for S Fig 1 panel C
├── 12_make_results_table.py           collates artifacts into Tables 2-4
│
├── figures/
│   ├── fig_style.py                   shared matplotlib style + helpers
│   ├── 13_make_fig1.py                Figure 1: study design schematic
│   ├── 14_make_fig2.py                Figure 2: locked 7-gene panel
│   ├── 15_make_fig3.py                Figure 3: internal validation on GPL570
│   ├── 16_make_fig4.py                Figure 4: GPL6244 external validation
│   ├── 17_make_fig5.py                Figure 5: GSE42133 + GSE25507 ROC
│   ├── 18_make_fig6.py                Figure 6: monocyte tracking
│   ├── 19_make_graphical_abstract.py  graphical abstract (SVG + PNG + PDF)
│   ├── 20_make_sfig1.py               S Fig 1: 7-gene annotation landscape
│   └── 21_make_sfig2.py               S Fig 2: LASSO regularisation path
│
└── tables/
    ├── 22_make_main_tables.py         Tables 1-4 (three-line, black-and-white)
    └── 23_make_supplementary_tables.py Tables S1-S6 (three-line, black-and-white)
```

## Output deliverables produced by this code

| Deliverable | Producing script |
|---|---|
| `figures/output/Figure_1.{png,pdf,tiff}` | `13_make_fig1.py` |
| `figures/output/Figure_2.{png,pdf,tiff}` | `14_make_fig2.py` |
| `figures/output/Figure_3.{png,pdf,tiff}` | `15_make_fig3.py` |
| `figures/output/Figure_4.{png,pdf,tiff}` | `16_make_fig4.py` |
| `figures/output/Figure_5.{png,pdf,tiff}` | `17_make_fig5.py` |
| `figures/output/Figure_6.{png,pdf,tiff}` | `18_make_fig6.py` |
| `figures/output/Graphical_Abstract.{png,pdf,svg}` | `19_make_graphical_abstract.py` |
| `figures/output/Supplementary_Figure_1.{png,pdf,tiff}` | `20_make_sfig1.py` |
| `figures/output/Supplementary_Figure_2.{png,pdf,tiff}` | `21_make_sfig2.py` |
| `Tables_revised.docx` (4 main tables) | `tables/22_make_main_tables.py` |
| `Supplementary_Tables.docx` (6 supp tables) | `tables/23_make_supplementary_tables.py` |

All TIFFs at 300 dpi, LZW compression.  PDFs with embedded TrueType
fonts (fonttype 42).  Tables in black-and-white three-line format with
Times New Roman.

## Sanity-check (numbers this code must reproduce)

| Cohort                        | Code | Manuscript |
|---|---|---|
| GSE18123-GPL570 (apparent)    | 0.9421 | 0.9421 |
| Nested CV mean                | 0.762  | 0.762  |
| Bootstrap .632 corrected      | 0.891  | 0.891  |
| GSE18123-GPL6244 adapted      | 0.6890 | 0.6890 |
| GPL6244 males                 | 0.7044 | 0.7047 |
| GPL6244 females               | 0.6397 | 0.6342 |
| GSE42133 6-gene refit         | 0.4195 | 0.4195 |
| GSE42133 7-gene zero-fill     | 0.4209 | 0.4209 |
| GSE25507 adapted              | 0.4636 | 0.4636 |
| GSE25507 within-tissue refit  | 0.620  | 0.620  |
| Monocyte Pearson r (4 cohorts) | +0.41 to +0.57 |

All numbers match to four decimal places.  Random seeds (= 42) are set
in `config.py`.

## Software versions tested

- Python 3.10 / 3.12
- numpy ≥ 1.24, pandas ≥ 2.0, scipy ≥ 1.10
- scikit-learn ≥ 1.3, statsmodels ≥ 0.14
- matplotlib ≥ 3.7, pyarrow ≥ 14, python-docx ≥ 1.1

## License

MIT (see `LICENSE`).

## Citation

> Wen J. et al. (2026).  A monocyte-tracking blood expression score reveals
> platform-specific dissociation between cellular composition and autism
> case-control status.  *iScience*, in revision.
>
> Code archive: Zenodo DOI 10.5281/zenodo.20404078
