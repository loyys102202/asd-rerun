#!/usr/bin/env bash
# ===========================================================================
# run_all.sh — end-to-end reproduction of ISCIENCE-D-26-03472 deliverables.
#
# Two modes:
#   bash run_all.sh         → reproduce from raw GEO (downloads ~1 GB)
#   bash run_all.sh figures → assume artifacts present, just regenerate
#                             figures + tables (fast, ~5 min)
#
# Random seed = 42 throughout (see config.py).
# ===========================================================================
set -e
MODE=${1:-full}

if [ "$MODE" = "full" ]; then
  echo "=== [Stage 1] Download GEO data ==="
  python3 01_download_geo.py

  echo "=== [Stage 2] Preprocess (log2 + normalise) ==="
  python3 02_preprocess.py

  echo "=== [Stage 3] Common-gene pool ==="
  python3 03_common_genes.py

  echo "=== [Stage 4] Nested CV with stability selection ==="
  python3 04_nested_cv.py

  echo "=== [Stage 5] Lock 7-gene panel ==="
  python3 05_lock_panel.py

  echo "=== [Stage 6] Bootstrap optimism (.632) ==="
  python3 06_bootstrap.py

  echo "=== [Stage 7] Training permutation null ==="
  python3 07_permutation_train.py

  echo "=== [Stage 8] External validation: GSE18123-GPL6244 ==="
  python3 08_external_GPL6244.py

  echo "=== [Stage 9] Independent external: GSE42133 ==="
  python3 09_external_GSE42133.py

  echo "=== [Stage 10] Tissue boundary: GSE25507 ==="
  python3 10_tissue_boundary_GSE25507.py

  echo "=== [Stage 11] Cell-type marker correlations ==="
  python3 11_cell_markers.py

  echo "=== [Stage 11b] Bootstrap fold change for S Fig 1 ==="
  python3 11b_precompute_fold_change.py

  echo "=== [Stage 12] Collate results table ==="
  python3 12_make_results_table.py
fi

# Figures (always rerun in both modes)
echo ""
echo "=== [Stage 13-21] Generate all figures ==="
cd figures
python3 13_make_fig1.py
python3 14_make_fig2.py
python3 15_make_fig3.py
python3 16_make_fig4.py
python3 17_make_fig5.py
python3 18_make_fig6.py
python3 19_make_graphical_abstract.py
python3 20_make_sfig1.py
python3 21_make_sfig2.py
cd ..

# Tables (always rerun in both modes)
echo ""
echo "=== [Stage 22-23] Generate tables ==="
cd tables
python3 22_make_main_tables.py
python3 23_make_supplementary_tables.py
cd ..

echo ""
echo "All outputs ready:"
echo "  Figures: figures/output/Figure_{1..6}.{png,pdf,tiff}"
echo "           figures/output/Supplementary_Figure_{1,2}.{png,pdf,tiff}"
echo "           figures/output/Graphical_Abstract.{png,pdf,svg}"
echo "  Tables:  Tables_revised.docx"
echo "           Supplementary_Tables.docx"
