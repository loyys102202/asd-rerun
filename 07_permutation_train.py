"""07_permutation_train.py
================================================================================
Label-permutation null testing on the GSE18123-GPL570 training cohort.

For each of 1,000 random permutations of the ASD / control labels:
  1. Refit the LASSO LP on the locked panel under the permuted labels.
  2. Compute the apparent (training-fit) AUC.
The empirical permutation p-value is (1 + sum(perm AUC ≥ observed)) / (n_perms + 1).

This procedure holds the gene panel fixed and asks how often the apparent
training-fit AUC under a randomly relabelled training set would meet or
exceed the observed apparent AUC on the true labelling.

Inputs
------
data/processed/expr_gene_GPL570.parquet
data/processed/phenotype_GPL570.csv
data/artifacts/locked_model_main_7gene.json

Outputs
-------
appends 'permutation_training' section to data/artifacts/robustness_checks.json
"""
from __future__ import annotations
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

import config
from pipeline import read_json, write_json, get_lasso_cs


def main():
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GPL570.csv")
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GPL570.parquet")
    model = read_json(config.LOCKED_MODEL_MAIN)

    X_panel = expr.loc[pheno["sample_id"].values, model["genes"]]
    y_true = pheno.set_index("sample_id").loc[X_panel.index, "asd"].astype(int).values

    # Standardise once (mean/SD from full training set)
    mu = X_panel.mean()
    sd = X_panel.std()
    sd[sd == 0] = 1.0
    Xz = ((X_panel - mu) / sd).values

    # Observed apparent AUC under true labels
    cs = get_lasso_cs()
    inner = StratifiedKFold(n_splits=config.INNER_CV_FOLDS, shuffle=True,
                            random_state=config.SEED)
    clf = LogisticRegressionCV(Cs=cs, penalty="l1", solver="liblinear",
                                cv=inner, scoring="roc_auc",
                                max_iter=config.LASSO_MAX_ITER,
                                random_state=config.SEED)
    clf.fit(Xz, y_true)
    observed = float(roc_auc_score(y_true, clf.predict_proba(Xz)[:, 1]))
    print(f"observed apparent AUC: {observed:.4f}")

    rng = np.random.RandomState(config.SEED)
    n_perms = config.N_PERMS_TRAIN
    n_ge = 0
    null_aucs = []
    y_perm = y_true.copy()
    for k in range(n_perms):
        rng.shuffle(y_perm)
        cv_p = StratifiedKFold(n_splits=config.INNER_CV_FOLDS, shuffle=True,
                               random_state=int(rng.randint(2**31 - 1)))
        try:
            clf_p = LogisticRegressionCV(Cs=cs, penalty="l1", solver="liblinear",
                                          cv=cv_p, scoring="roc_auc",
                                          max_iter=config.LASSO_MAX_ITER)
            clf_p.fit(Xz, y_perm)
            auc_p = float(roc_auc_score(y_perm, clf_p.predict_proba(Xz)[:, 1]))
        except Exception:
            continue
        null_aucs.append(auc_p)
        if auc_p >= observed:
            n_ge += 1
        if (k + 1) % 100 == 0:
            print(f"  perm {k + 1}/{n_perms} done")

    p_value = (n_ge + 1) / (len(null_aucs) + 1)
    result = {
        "observed_apparent_auc": observed,
        "n_permutations": len(null_aucs),
        "n_ge_observed": n_ge,
        "p_value": float(p_value),
        "null_mean": float(np.mean(null_aucs)),
        "null_sd": float(np.std(null_aucs, ddof=1)),
        "null_p95": float(np.percentile(null_aucs, 95)),
    }

    print(f"\npermutation null AUC mean ± SD: {result['null_mean']:.4f} ± {result['null_sd']:.4f}")
    print(f"empirical p-value: {result['p_value']:.4f}")

    if config.ROBUSTNESS_CHECKS.exists():
        store = read_json(config.ROBUSTNESS_CHECKS)
    else:
        store = {}
    store["permutation_training"] = result
    write_json(store, config.ROBUSTNESS_CHECKS)
    print(f"\nwrote {config.ROBUSTNESS_CHECKS}")
    print("stage 07 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
