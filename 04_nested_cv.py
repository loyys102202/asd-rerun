"""04_nested_cv.py
================================================================================
Repeated nested cross-validation with t-test pre-screen (top 500 genes) and
LASSO logistic regression, on the GSE18123-GPL570 training cohort.

Design
------
- Outer: 5-fold stratified CV, repeated 100 times → 500 outer folds.
- Within each outer fold:
    * Welch's two-sample t-test on training subset → keep top 500 genes
    * Standardise top-500 features using training-subset mean and SD
    * Inner 5-fold stratified CV LASSO via LogisticRegressionCV
      over Cs = logspace(-3, 1, 30)
    * Predict on held-out outer fold; record fold AUC and selected genes
- Selection frequency per gene = (folds with non-zero coef) / 500.
- Per-repeat mean AUC = mean of the 5 outer-fold AUCs in that repeat;
  primary internal estimator = mean of 100 per-repeat means,
  95% CI = empirical 2.5–97.5 percentiles of those 100 means.

Inputs
------
data/processed/expr_gene_GPL570.parquet
data/processed/phenotype_GPL570.csv
data/processed/common_genes_570_6244.txt

Outputs
-------
data/artifacts/nestedcv_outer_folds.csv  (per-fold: repeat, fold, auc, n_selected, selected_genes)
data/artifacts/nestedcv_gene_freq.csv    (gene, count, frequency) sorted descending
data/artifacts/nestedcv_auc_by_rep.csv   (repeat, mean_auc)
"""
from __future__ import annotations
import sys
import json
from collections import Counter

import numpy as np
import pandas as pd

import config
from pipeline import nested_cv_lasso_stability, get_lasso_cs


def main(n_jobs: int = -1):
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GPL570.parquet")
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GPL570.csv")
    common = config.COMMON_GENES.read_text().strip().split("\n")

    # Align samples & restrict to common gene pool
    common_in_expr = [g for g in common if g in expr.columns]
    if len(common_in_expr) != len(common):
        print(f"NOTE: {len(common) - len(common_in_expr)} common genes absent from GPL570 expression "
              f"(probably probe→gene aggregation differences); using {len(common_in_expr)}")
    X = expr.loc[pheno["sample_id"].values, common_in_expr]
    y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int)

    print(f"training cohort: n_samples={len(X)}, n_features={X.shape[1]}, "
          f"n_ASD={int(y.sum())}, n_ctrl={int((1 - y).sum())}")

    lasso_cs = get_lasso_cs()
    print(f"LASSO grid: {len(lasso_cs)} values from {lasso_cs[0]:.2e} to {lasso_cs[-1]:.2e}")

    results = nested_cv_lasso_stability(
        X, y,
        n_outer=config.N_OUTER_FOLDS,
        n_repeats=config.N_REPEATS,
        top_k=config.TOP_K_TTEST,
        inner_folds=config.INNER_CV_FOLDS,
        lasso_cs=lasso_cs,
        lasso_max_iter=config.LASSO_MAX_ITER,
        seed=config.SEED,
        n_jobs=n_jobs,
    )

    # Per-fold table
    rows = []
    for r in results:
        rows.append({
            "repeat": r["repeat"],
            "fold": r["fold"],
            "auc": r["auc"],
            "n_selected": r["n_selected"],
            "selected_genes": ";".join(r["selected_genes"]),
        })
    fold_df = pd.DataFrame(rows)
    fold_df.to_csv(config.NESTEDCV_FOLD_RESULTS, index=False)
    print(f"\nwrote {config.NESTEDCV_FOLD_RESULTS}")

    # Per-repeat mean AUC
    per_rep = fold_df.groupby("repeat")["auc"].mean().reset_index()
    per_rep.columns = ["repeat", "mean_auc"]
    per_rep.to_csv(config.NESTEDCV_AUC_BY_REP, index=False)

    rep_means = per_rep["mean_auc"].values
    auc_mean = float(rep_means.mean())
    auc_low = float(np.percentile(rep_means, 2.5))
    auc_high = float(np.percentile(rep_means, 97.5))
    print(f"nested CV AUC: {auc_mean:.4f}  95% CI [{auc_low:.4f}, {auc_high:.4f}]")

    # Gene selection frequencies
    n_folds = len(results)
    counts = Counter()
    for r in results:
        for g in r["selected_genes"]:
            counts[g] += 1
    freq_df = pd.DataFrame([
        {"gene": g, "count": c, "frequency": c / n_folds}
        for g, c in counts.items()
    ]).sort_values("frequency", ascending=False).reset_index(drop=True)
    freq_df.to_csv(config.NESTEDCV_GENE_FREQ, index=False)
    print(f"wrote {config.NESTEDCV_GENE_FREQ}  ({len(freq_df)} genes ever selected)")

    print(f"\nGenes at or above {int(config.STABILITY_THRESHOLD * 100)}% selection frequency:")
    print(freq_df[freq_df["frequency"] >= config.STABILITY_THRESHOLD].to_string(index=False))

    print("\nstage 04 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
