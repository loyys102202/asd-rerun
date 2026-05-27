"""06_bootstrap.py
================================================================================
0.632 bootstrap optimism correction for the locked 7-gene panel.

Procedure
---------
For B = 1000 bootstrap resamples of the training cohort:
  1. Refit the LASSO LP on the locked panel using the bootstrap sample
     (standardised by the bootstrap sample's own mean/SD).
  2. Compute apparent AUC on the bootstrap sample.
  3. Compute test AUC on the original full sample.
  4. Optimism = apparent_b − test_b.
Corrected AUC = apparent_full − mean(optimism).

Inputs
------
data/processed/expr_gene_GPL570.parquet
data/processed/phenotype_GPL570.csv
data/artifacts/locked_model_main_7gene.json

Outputs
-------
appends 'bootstrap' section to data/artifacts/robustness_checks.json
"""
from __future__ import annotations
import sys

import pandas as pd

import config
from pipeline import (bootstrap_optimism, read_json, write_json, get_lasso_cs)


def main():
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GPL570.csv")
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GPL570.parquet")
    model = read_json(config.LOCKED_MODEL_MAIN)

    X = expr.loc[pheno["sample_id"].values, model["genes"]]
    y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int).values

    print(f"bootstrap optimism: B = {config.B_BOOTSTRAP}, n = {len(y)}")
    result = bootstrap_optimism(
        X_train=X, y_train=y,
        panel_genes=model["genes"],
        lasso_cs=get_lasso_cs(),
        inner_folds=config.INNER_CV_FOLDS,
        lasso_max_iter=config.LASSO_MAX_ITER,
        B=config.B_BOOTSTRAP,
        seed=config.SEED,
    )

    print(f"  apparent AUC        : {result['apparent_auc']:.4f}")
    print(f"  mean optimism       : {result['mean_optimism']:.4f}")
    print(f"  corrected AUC       : {result['optimism_corrected_auc']:.4f}")
    print(f"  B used / requested  : {result['B_used']} / {result['B_requested']}")

    # Merge into robustness_checks.json
    if config.ROBUSTNESS_CHECKS.exists():
        store = read_json(config.ROBUSTNESS_CHECKS)
    else:
        store = {}
    store["bootstrap_optimism"] = result
    write_json(store, config.ROBUSTNESS_CHECKS)
    print(f"\nwrote {config.ROBUSTNESS_CHECKS}")
    print("stage 06 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
