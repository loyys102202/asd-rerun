"""11_cell_markers.py
================================================================================
Per-cohort Pearson correlations of LP7 with marker-score signatures for
eight blood cell types. The marker lists are compiled from CIBERSORT LM22
(Newman et al. 2015) and xCell (Aran et al. 2017) and are defined in
pipeline.CELL_MARKERS.

For each cohort and each cell type:
  1. Compute per-sample marker score = mean of within-cohort z-scored
     expression across that cell type's marker genes.
  2. Compute Pearson r between marker score and LP7.

Cell-composition-adjusted external AUC on GPL6244 is also reported:
  for each cell type, regress LP7 on the marker score, then recompute
  AUC using residuals (this tests whether ASD discrimination survives
  removal of that cell type's signal).

Inputs
------
data/processed/expr_gene_<KEY>.parquet for KEY ∈ {GPL570, GPL6244, GSE25507, GSE42133}
data/processed/phenotype_<KEY>.csv
data/artifacts/locked_model_main_7gene.json

Outputs
-------
data/artifacts/cell_marker_analysis.json
data/artifacts/cell_marker_GSE42133.json (the GSE42133 portion, for cross-check)
"""
from __future__ import annotations
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

import config
from pipeline import (read_json, write_json, predict_proba_locked,
                       cell_marker_correlations, CELL_MARKERS,
                       cell_type_marker_score)


def _compute_lp7(model: dict, X: pd.DataFrame, mode: str) -> pd.Series:
    """Compute LP7 (linear predictor on logit scale) using the given mode."""
    # We need the LP itself (logit), not the probability, because the LP is
    # what gets correlated with cell-marker scores in the paper.
    coefs = pd.Series(model["coef"], index=model["genes"])
    intercept = float(model["intercept"])
    present = [g for g in model["genes"] if g in X.columns]

    if mode == "platform_adapted":
        Z = (X[present] - X[present].mean()) / X[present].std()
    elif mode == "strict_train":
        mu = pd.Series(model["train_mean"], index=model["genes"])
        sd = pd.Series(model["train_sd"], index=model["genes"])
        Z = (X[present] - mu[present]) / sd[present]
    else:
        raise ValueError(f"unknown mode {mode!r}")

    lp = intercept + (Z * coefs[present]).sum(axis=1)
    return lp


def main():
    model = read_json(config.LOCKED_MODEL_MAIN)

    out = {}
    out_gse42133 = {}

    for key in config.ALL_KEYS:
        expr = pd.read_parquet(config.PROCESSED / f"expr_gene_{key}.parquet")
        pheno = pd.read_csv(config.PROCESSED / f"phenotype_{key}.csv")
        sample_ids = pheno["sample_id"].values
        X = expr.loc[sample_ids]
        y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int).values

        # Compute LP7 (use platform-adapted for all non-training cohorts;
        # for training, use the training-set standardisation directly so the
        # correlation reflects the model on its own scale)
        if key == "GPL570":
            lp7 = _compute_lp7(model, X, mode="strict_train")
        else:
            lp7 = _compute_lp7(model, X, mode="platform_adapted")

        cohort_result = cell_marker_correlations(X, lp7)
        n_ASD = int((y == 1).sum())
        n_CTL = int((y == 0).sum())
        out[key] = {
            "n_total": len(y), "n_ASD": n_ASD, "n_ctrl": n_CTL,
            "lp7_mean": float(lp7.mean()),
            "lp7_sd": float(lp7.std()),
            "correlations": cohort_result,
        }
        if key == "GSE42133":
            out_gse42133 = {key: out[key]}

        print(f"\n[{key}]  n={len(y)}  ASD={n_ASD}  ctrl={n_CTL}")
        for ct, info in cohort_result.items():
            print(f"  {ct:18s}  r = {info['r']:+.4f}  p = {info['p']:.2e}  "
                  f"(n_markers = {info['n_markers']})")

    # Cell-composition-adjusted external AUC on GPL6244
    print("\ncell-composition-adjusted external AUC on GPL6244 (per cell type):")
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GPL6244.parquet")
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GPL6244.csv")
    X = expr.loc[pheno["sample_id"].values]
    y = pheno.set_index("sample_id").loc[X.index, "asd"].astype(int).values
    p_adapt = predict_proba_locked(model, X, mode="platform_adapted").values
    base_auc = float(roc_auc_score(y, p_adapt))
    print(f"  baseline AUC: {base_auc:.4f}")

    adj_aucs = {}
    for ct, markers in CELL_MARKERS.items():
        score = cell_type_marker_score(X, markers).values
        if np.isnan(score).all():
            adj_aucs[ct] = {"adjusted_auc": None, "delta": None}
            continue
        # Residualise LP (use log-odds to avoid sigmoid compression)
        eps = 1e-9
        lp_adapt = np.log(np.clip(p_adapt, eps, 1 - eps) /
                          (1 - np.clip(p_adapt, eps, 1 - eps)))
        valid = ~np.isnan(score)
        reg = LinearRegression()
        reg.fit(score[valid].reshape(-1, 1), lp_adapt[valid])
        resid = lp_adapt.copy()
        resid[valid] = lp_adapt[valid] - reg.predict(score[valid].reshape(-1, 1))
        a = float(roc_auc_score(y, resid))
        adj_aucs[ct] = {"adjusted_auc": a, "delta": a - base_auc}
        print(f"  adjust for {ct:18s}: AUC {a:.4f}  ΔAUC {a - base_auc:+.4f}")

    out["cell_composition_adjusted_GPL6244"] = {
        "baseline_auc": base_auc,
        "adjusted": adj_aucs,
    }

    write_json(out, config.CELL_MARKER_ANALYSIS)
    write_json(out_gse42133, config.CELL_MARKER_GSE42133)
    print(f"\nwrote {config.CELL_MARKER_ANALYSIS}")
    print(f"wrote {config.CELL_MARKER_GSE42133}")
    print("stage 11 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
