"""Paths, seeds, and hyperparameters shared across all stages.

Importing this module also creates the data directories if they do not exist.
"""
from pathlib import Path

# ------------------------------------------------------------------------
# Directory layout
# ------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RAW = DATA / "raw"           # downloaded GEO series matrices and platform tables
PROCESSED = DATA / "processed"  # gene-aggregated expression matrices and phenotypes
ARTIFACTS = DATA / "artifacts"  # results: locked model, AUCs, calibration, etc.

for _d in (DATA, RAW, PROCESSED, ARTIFACTS):
    _d.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------
# Cohorts and their roles in the analysis
# ------------------------------------------------------------------------
# key -> {gse, platform, tissue, role}
COHORTS = {
    "GPL570": {
        "gse": "GSE18123",
        "platform": "GPL570",
        "tissue": "whole_blood",
        "role": "training",
    },
    "GPL6244": {
        "gse": "GSE18123",
        "platform": "GPL6244",
        "tissue": "whole_blood",
        "role": "external_primary",
    },
    "GSE25507": {
        "gse": "GSE25507",
        "platform": "GPL570",
        "tissue": "lymphocyte_enriched",
        "role": "tissue_boundary",
    },
    "GSE42133": {
        "gse": "GSE42133",
        "platform": "GPL10558",
        "tissue": "leukocyte",
        "role": "external_independent",
    },
}

# Convenience views
TRAINING_KEY = "GPL570"
EXTERNAL_KEYS = ["GPL6244", "GSE42133"]
TISSUE_BOUNDARY_KEY = "GSE25507"
ALL_KEYS = list(COHORTS.keys())

# ------------------------------------------------------------------------
# Pipeline hyperparameters
# ------------------------------------------------------------------------
# Master random seed — every stochastic procedure derives its RNG from this
SEED = 42

# Nested cross-validation
N_OUTER_FOLDS = 5
N_REPEATS = 100              # outer-CV repeats → 500 outer folds total
TOP_K_TTEST = 500            # genes retained by the within-fold t-test pre-screen
INNER_CV_FOLDS = 5
LASSO_LOG10_C_MIN = -3       # Cs = np.logspace(-3, 1, 30)
LASSO_LOG10_C_MAX = 1
LASSO_N_CS = 30
LASSO_MAX_ITER = 4000
STABILITY_THRESHOLD = 0.60   # gene kept in panel if selected in >= 60% of outer folds

# Bootstrap optimism correction
B_BOOTSTRAP = 1000

# Label-permutation null testing
N_PERMS_TRAIN = 1000
N_PERMS_EXT = 10000

# Locked panel constants (computed by 05_lock_panel.py; recorded here for
# downstream scripts to cross-check that the lock has not drifted unexpectedly)
EXPECTED_PANEL_GENES = ["KYNU", "CMYA5", "MAPK8IP1", "CES1",
                        "KIAA0087", "POU2AF1", "ERVK3-2"]

# Genes absent on GPL10558 (Illumina HT-12 v4) — used by 09_external_GSE42133.py
GPL10558_MISSING_PANEL_GENES = ["ERVK3-2"]
REDUCED_PANEL_FOR_GSE42133 = [g for g in EXPECTED_PANEL_GENES
                              if g not in GPL10558_MISSING_PANEL_GENES]

# ------------------------------------------------------------------------
# File paths for key artifacts
# ------------------------------------------------------------------------
LOCKED_MODEL_MAIN = ARTIFACTS / "locked_model_main_7gene.json"
LOCKED_MODEL_REDUCED = ARTIFACTS / "locked_model_reduced_6gene.json"
NESTEDCV_FOLD_RESULTS = ARTIFACTS / "nestedcv_outer_folds.csv"
NESTEDCV_GENE_FREQ = ARTIFACTS / "nestedcv_gene_freq.csv"
NESTEDCV_AUC_BY_REP = ARTIFACTS / "nestedcv_auc_by_rep.csv"
ROBUSTNESS_CHECKS = ARTIFACTS / "robustness_checks.json"
EXTERNAL_GPL6244 = ARTIFACTS / "external_validation_GPL6244.csv"
EXTERNAL_GSE42133 = ARTIFACTS / "external_validation_GSE42133.json"
GSE25507_EXPLORATION = ARTIFACTS / "GSE25507_exploration.json"
CALIBRATION = ARTIFACTS / "calibration.json"
CELL_MARKER_ANALYSIS = ARTIFACTS / "cell_marker_analysis.json"
CELL_MARKER_GSE42133 = ARTIFACTS / "cell_marker_GSE42133.json"
RESULTS_SUMMARY = ARTIFACTS / "results_summary.json"
COMMON_GENES = PROCESSED / "common_genes_570_6244.txt"
