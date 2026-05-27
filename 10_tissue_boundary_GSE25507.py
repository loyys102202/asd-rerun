"""10_tissue_boundary_GSE25507.py
================================================================================
Exploratory tissue-boundary analysis on GSE25507 (lymphocyte-enriched
mononuclear cells; n = 146). This cohort is NOT used as external validation
because its blood-fraction preparation differs substantially from the
training tissue (whole blood).

Reports
-------
- Platform-adapted standardisation: AUC, CI, two-sided permutation p
- Strict training-scale standardisation: AUC (sensitivity)
- Within-tissue refit: holding panel genes fixed, refit coefficients on the
  GSE25507 cohort itself, then report apparent AUC and direction concordance
  against the whole-blood-trained coefficients

No sign-flipping is performed. The original-direction AUC is reported as-is.

Inputs
------
data/processed/expr_gene_GSE25507.parquet
data/processed/phenotype_GSE25507.csv
data/artifacts/locked_model_main_7gene.json

Outputs
-------
data/artifacts/GSE25507_exploration.json
"""
from __future__ import annotations
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

import config
from pipeline import (read_json, write_json, predict_proba_locked,
                       auc as compute_auc, delong_ci, permutation_pvalue,
                       get_lasso_cs)


def main():
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GSE25507.csv")
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GSE25507.parquet")
    model = read_json(config.LOCKED_MODEL_MAIN)

    X = expr.loc[pheno["sample_id"].values]
    y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int).values

    panel_present = [g for g in model["genes"] if g in X.columns]
    panel_missing = [g for g in model["genes"] if g not in X.columns]
    n_ASD = int((y == 1).sum())
    n_CTL = int((y == 0).sum())
    print(f"GSE25507: n_total = {len(y)}, n_ASD = {n_ASD}, n_ctrl = {n_CTL}")
    print(f"  panel genes present: {panel_present}")
    if panel_missing:
        print(f"  panel genes MISSING: {panel_missing}")

    out = {"n_total": len(y), "n_ASD": n_ASD, "n_ctrl": n_CTL,
           "panel_present": panel_present, "panel_missing": panel_missing}

    # ---------- Platform-adapted (primary) ----------
    p_adapt = predict_proba_locked(model, X, mode="platform_adapted").values
    auc_adapt = compute_auc(y, p_adapt)
    ci_adapt = delong_ci(y, p_adapt)
    perm_p = permutation_pvalue(y, p_adapt, n_perms=config.N_PERMS_EXT,
                                  seed=config.SEED, two_sided=True)
    print(f"\nplatform-adapted AUC: {auc_adapt:.4f}  CI [{ci_adapt[0]:.4f}, {ci_adapt[1]:.4f}]")
    print(f"two-sided perm p:      {perm_p:.4f}")
    out["platform_adapted"] = {
        "auc": auc_adapt, "ci_low": ci_adapt[0], "ci_high": ci_adapt[1],
        "perm_p_two_sided": perm_p,
    }

    # ---------- Strict training-scale (sensitivity) ----------
    p_strict = predict_proba_locked(model, X, mode="strict_train").values
    auc_strict = compute_auc(y, p_strict)
    print(f"strict train-scale AUC: {auc_strict:.4f}")
    out["strict_train_scale"] = {"auc": auc_strict}

    # ---------- Within-tissue refit ----------
    # Hold panel genes fixed, refit coefficients on this cohort.
    print("\nwithin-tissue refit on GSE25507...")
    X_panel = X[panel_present]
    mu = X_panel.mean()
    sd = X_panel.std()
    sd[sd == 0] = 1.0
    Xz = ((X_panel - mu) / sd).values

    cs = get_lasso_cs()
    inner = StratifiedKFold(n_splits=config.INNER_CV_FOLDS, shuffle=True,
                            random_state=config.SEED)
    clf = LogisticRegressionCV(Cs=cs, penalty="l1", solver="liblinear",
                                cv=inner, scoring="roc_auc",
                                max_iter=config.LASSO_MAX_ITER,
                                random_state=config.SEED)
    clf.fit(Xz, y)
    refit_auc = float(roc_auc_score(y, clf.predict_proba(Xz)[:, 1]))
    refit_coefs = dict(zip(panel_present, [float(c) for c in clf.coef_[0]]))
    refit_intercept = float(clf.intercept_[0])

    # Direction concordance
    train_coefs = dict(zip(model["genes"], model["coef"]))
    concord = 0
    direction_table = []
    for g in panel_present:
        s_train = np.sign(train_coefs[g])
        s_refit = np.sign(refit_coefs[g])
        same = bool(s_train == s_refit and s_train != 0 and s_refit != 0)
        if same:
            concord += 1
        direction_table.append({"gene": g,
                                  "train_coef": train_coefs[g],
                                  "refit_coef": refit_coefs[g],
                                  "concordant_direction": same})
    print(f"refit apparent AUC: {refit_auc:.4f}")
    print(f"direction concordance: {concord}/{len(panel_present)}")
    out["within_tissue_refit"] = {
        "panel_used": panel_present,
        "refit_apparent_auc": refit_auc,
        "refit_intercept": refit_intercept,
        "refit_coefs": refit_coefs,
        "direction_concordance_n": concord,
        "direction_concordance_d": len(panel_present),
        "direction_table": direction_table,
    }

    write_json(out, config.GSE25507_EXPLORATION)
    print(f"\nwrote {config.GSE25507_EXPLORATION}")
    print("stage 10 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
