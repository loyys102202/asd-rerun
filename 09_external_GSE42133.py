"""09_external_GSE42133.py
================================================================================
Independent external test on GSE42133 (Pramparo et al. 2015; Illumina HT-12 v4
toddler leukocyte; n = 147, all male).

ERVK3-2 has no reliable probe on Illumina HT-12 v4 (GPL10558). Two
complementary configurations are evaluated:
  - reduced 6-gene model (refit by 05_lock_panel.py and stored in
    locked_model_reduced_6gene.json)
  - 7-gene zero-fill sensitivity (substitute z = 0 for the missing gene
    while holding all other coefficients fixed)

Reports
-------
For each configuration: AUC, 95% CI (DeLong), Brier, two-sided permutation p.

Inputs
------
data/processed/expr_gene_GSE42133.parquet
data/processed/phenotype_GSE42133.csv
data/artifacts/locked_model_main_7gene.json
data/artifacts/locked_model_reduced_6gene.json

Outputs
-------
data/artifacts/external_validation_GSE42133.json
"""
from __future__ import annotations
import sys

import numpy as np
import pandas as pd

import config
from pipeline import (read_json, write_json, predict_proba_locked,
                       auc as compute_auc, brier as compute_brier,
                       delong_ci, permutation_pvalue)


def main():
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GSE42133.csv")
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GSE42133.parquet")
    main_model = read_json(config.LOCKED_MODEL_MAIN)
    reduced_model = read_json(config.LOCKED_MODEL_REDUCED)

    sample_ids = pheno["sample_id"].values
    X = expr.loc[sample_ids]
    y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int).values

    n_ASD = int((y == 1).sum())
    n_CTL = int((y == 0).sum())
    print(f"GSE42133: n_total = {len(y)}, n_ASD = {n_ASD}, n_ctrl = {n_CTL}")

    # Diagnostic: which panel genes are present?
    main_present = [g for g in main_model["genes"] if g in X.columns]
    main_missing = [g for g in main_model["genes"] if g not in X.columns]
    reduced_present = [g for g in reduced_model["genes"] if g in X.columns]
    reduced_missing = [g for g in reduced_model["genes"] if g not in X.columns]
    print(f"  7-gene panel present on GPL10558: {main_present}")
    print(f"  7-gene panel MISSING on GPL10558: {main_missing}")
    print(f"  6-gene reduced panel present:    {reduced_present}")
    print(f"  6-gene reduced panel MISSING:    {reduced_missing}")

    out = {"n_total": len(y), "n_ASD": n_ASD, "n_ctrl": n_CTL,
           "main_panel_present": main_present, "main_panel_missing": main_missing,
           "reduced_panel_present": reduced_present, "reduced_panel_missing": reduced_missing,
           "configurations": {}}

    # ---------- 6-gene refit (primary configuration) ----------
    print("\nconfiguration: 6-gene refit (ERVK3-2 dropped from linear predictor)")
    p_6gene = predict_proba_locked(reduced_model, X, mode="platform_adapted").values
    auc_6 = compute_auc(y, p_6gene)
    ci_6 = delong_ci(y, p_6gene)
    brier_6 = compute_brier(y, p_6gene)
    perm_p_6 = permutation_pvalue(y, p_6gene, n_perms=config.N_PERMS_EXT,
                                    seed=config.SEED, two_sided=True)
    print(f"  AUC: {auc_6:.4f}  CI [{ci_6[0]:.4f}, {ci_6[1]:.4f}]")
    print(f"  Brier: {brier_6:.4f}")
    print(f"  perm p (two-sided): {perm_p_6:.4f}")
    out["configurations"]["reduced_6gene"] = {
        "auc": auc_6, "ci_low": ci_6[0], "ci_high": ci_6[1],
        "brier": brier_6, "perm_p_two_sided": perm_p_6,
    }

    # ---------- 7-gene zero-fill sensitivity ----------
    print("\nconfiguration: 7-gene zero-fill (substitute z = 0 for ERVK3-2)")
    p_7zf = predict_proba_locked(main_model, X, mode="zero_fill").values
    auc_7 = compute_auc(y, p_7zf)
    ci_7 = delong_ci(y, p_7zf)
    brier_7 = compute_brier(y, p_7zf)
    perm_p_7 = permutation_pvalue(y, p_7zf, n_perms=config.N_PERMS_EXT,
                                    seed=config.SEED, two_sided=True)
    print(f"  AUC: {auc_7:.4f}  CI [{ci_7[0]:.4f}, {ci_7[1]:.4f}]")
    print(f"  Brier: {brier_7:.4f}")
    print(f"  perm p (two-sided): {perm_p_7:.4f}")
    out["configurations"]["zero_fill_7gene"] = {
        "auc": auc_7, "ci_low": ci_7[0], "ci_high": ci_7[1],
        "brier": brier_7, "perm_p_two_sided": perm_p_7,
    }

    write_json(out, config.EXTERNAL_GSE42133)
    print(f"\nwrote {config.EXTERNAL_GSE42133}")
    print("stage 09 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
