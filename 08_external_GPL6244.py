"""08_external_GPL6244.py
================================================================================
External validation on GSE18123-GPL6244 (n = 186) — cross-platform within
the same GEO series.

Reports
-------
- Platform-adapted standardisation: AUC, 95% CI (DeLong), Brier
- Strict training-scale standardisation: AUC, CI, Brier  (sensitivity)
- Sex-stratified AUC + CI for males and females (platform-adapted)
- Logistic recalibration: fit α, β; recalibrated Brier
- Label-permutation null (10,000 perms)

Inputs
------
data/processed/expr_gene_GPL6244.parquet
data/processed/phenotype_GPL6244.csv
data/artifacts/locked_model_main_7gene.json

Outputs
-------
data/artifacts/external_validation_GPL6244.csv  (long format, one row per analysis)
data/artifacts/calibration.json                 (α, β, recal Brier, etc.)
"""
from __future__ import annotations
import sys

import numpy as np
import pandas as pd

import config
from pipeline import (read_json, write_json, predict_proba_locked,
                       auc as compute_auc, brier as compute_brier,
                       delong_ci, permutation_pvalue, logistic_recalibration,
                       linear_predictor, sigmoid)


def main():
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GPL6244.csv")
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GPL6244.parquet")
    model = read_json(config.LOCKED_MODEL_MAIN)

    # Align samples
    X = expr.loc[pheno["sample_id"].values, [g for g in model["genes"] if g in expr.columns]]
    y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int).values

    panel = model["genes"]
    missing = [g for g in panel if g not in X.columns]
    if missing:
        print(f"WARNING: panel genes missing from GPL6244 expression: {missing}")
    n_ASD = int((y == 1).sum())
    n_CTL = int((y == 0).sum())
    print(f"GPL6244: n_total = {len(y)}, n_ASD = {n_ASD}, n_ctrl = {n_CTL}")

    rows = []

    # ---------- Platform-adapted standardisation (primary analysis) ----------
    p_adapted = predict_proba_locked(model, X, mode="platform_adapted").values
    lp_adapted = np.log(p_adapted / (1 - p_adapted))
    auc_adapt = compute_auc(y, p_adapted)
    ci_adapt = delong_ci(y, p_adapted)
    brier_adapt = compute_brier(y, p_adapted)
    rows.append({"analysis": "platform_adapted", "n_ASD": n_ASD, "n_ctrl": n_CTL,
                 "auc": auc_adapt, "ci_low": ci_adapt[0], "ci_high": ci_adapt[1],
                 "brier": brier_adapt})
    print(f"\nplatform-adapted   : AUC {auc_adapt:.4f}  CI [{ci_adapt[0]:.4f}, {ci_adapt[1]:.4f}]  Brier {brier_adapt:.4f}")

    # ---------- Strict training-scale standardisation (sensitivity) ----------
    p_strict = predict_proba_locked(model, X, mode="strict_train").values
    auc_strict = compute_auc(y, p_strict)
    ci_strict = delong_ci(y, p_strict)
    brier_strict = compute_brier(y, p_strict)
    rows.append({"analysis": "strict_train", "n_ASD": n_ASD, "n_ctrl": n_CTL,
                 "auc": auc_strict, "ci_low": ci_strict[0], "ci_high": ci_strict[1],
                 "brier": brier_strict})
    print(f"strict train scale : AUC {auc_strict:.4f}  CI [{ci_strict[0]:.4f}, {ci_strict[1]:.4f}]  Brier {brier_strict:.4f}")

    # ---------- Sex stratification (platform-adapted) ----------
    sex = pheno.set_index("sample_id").loc[X.index, "sex"]
    for sex_label, mask_label in [("M", "males"), ("F", "females")]:
        mask = (sex == sex_label).values
        if mask.sum() < 5 or len(np.unique(y[mask])) < 2:
            print(f"  {mask_label}: insufficient n for stratified analysis")
            continue
        p_sub = p_adapted[mask]
        y_sub = y[mask]
        a = compute_auc(y_sub, p_sub)
        ci = delong_ci(y_sub, p_sub)
        b = compute_brier(y_sub, p_sub)
        n_ASD_sub = int((y_sub == 1).sum())
        n_CTL_sub = int((y_sub == 0).sum())
        rows.append({"analysis": f"platform_adapted_{mask_label}",
                     "n_ASD": n_ASD_sub, "n_ctrl": n_CTL_sub,
                     "auc": a, "ci_low": ci[0], "ci_high": ci[1], "brier": b})
        print(f"{mask_label:9s}        : AUC {a:.4f}  CI [{ci[0]:.4f}, {ci[1]:.4f}]  "
              f"Brier {b:.4f}  (n_ASD={n_ASD_sub}, n_ctrl={n_CTL_sub})")

    # ---------- Permutation null on platform-adapted ----------
    print(f"\npermutation null ({config.N_PERMS_EXT} perms, two-sided)...")
    perm_p = permutation_pvalue(y, p_adapted, n_perms=config.N_PERMS_EXT,
                                 seed=config.SEED, two_sided=True)
    rows.append({"analysis": "platform_adapted_permutation_p", "n_ASD": n_ASD,
                 "n_ctrl": n_CTL, "auc": auc_adapt, "ci_low": np.nan, "ci_high": np.nan,
                 "brier": np.nan, "perm_p_two_sided": perm_p})
    print(f"two-sided perm p   : {perm_p:.6f}")

    df = pd.DataFrame(rows)
    df.to_csv(config.EXTERNAL_GPL6244, index=False)
    print(f"\nwrote {config.EXTERNAL_GPL6244}")

    # ---------- Logistic recalibration (platform-adapted LP) ----------
    # Use linear predictor on platform-adapted z-scores
    lp_for_recal = lp_adapted
    alpha, beta = logistic_recalibration(y, lp_for_recal)
    p_recal = sigmoid(alpha + beta * lp_for_recal)
    auc_recal = compute_auc(y, p_recal)
    brier_recal = compute_brier(y, p_recal)
    print(f"\nrecalibration: α = {alpha:+.4f},  β = {beta:+.4f}")
    print(f"  AUC unchanged: {auc_recal:.4f}")
    print(f"  Brier: {brier_adapt:.4f} → {brier_recal:.4f}")

    cal = {
        "platform_adapted": {
            "alpha": alpha, "beta": beta,
            "auc_pre": auc_adapt, "auc_post": auc_recal,
            "brier_pre": brier_adapt, "brier_post": brier_recal,
        }
    }
    write_json(cal, config.CALIBRATION)
    print(f"wrote {config.CALIBRATION}")
    print("\nstage 08 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
