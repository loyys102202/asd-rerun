"""05_lock_panel.py
================================================================================
Apply the ≥ 60% stability threshold to the nested-CV gene selection
frequencies, lock the resulting gene set, and refit final coefficients on
the full GSE18123-GPL570 training set.

Inputs
------
data/processed/expr_gene_GPL570.parquet
data/processed/phenotype_GPL570.csv
data/artifacts/nestedcv_gene_freq.csv

Outputs
-------
data/artifacts/locked_model_main_7gene.json
"""
from __future__ import annotations
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, brier_score_loss

import config
from pipeline import (refit_locked_panel, write_json, get_lasso_cs,
                       predict_proba_locked)


def main():
    freq_df = pd.read_csv(config.NESTEDCV_GENE_FREQ)
    panel = freq_df[freq_df["frequency"] >= config.STABILITY_THRESHOLD]["gene"].tolist()
    print(f"genes at or above {int(config.STABILITY_THRESHOLD * 100)}%: {len(panel)}")
    print(panel)

    if not panel:
        raise RuntimeError("empty panel — adjust stability threshold or check 04 output")

    # Remove collinear pseudogenes (e.g. CES1P1 ↔ CES1). Detect via pairwise r > 0.999.
    expr = pd.read_parquet(config.PROCESSED / "expr_gene_GPL570.parquet")
    pheno = pd.read_csv(config.PROCESSED / "phenotype_GPL570.csv")
    X_panel = expr.loc[pheno["sample_id"].values, [g for g in panel if g in expr.columns]]
    cmat = X_panel.corr().abs()
    # Find pairs with r > 0.999 (excluding self)
    dropped = []
    keep = list(X_panel.columns)
    for i, gi in enumerate(X_panel.columns):
        for gj in X_panel.columns[i + 1:]:
            if cmat.loc[gi, gj] > 0.999:
                # Drop the lower-frequency partner of the pair
                f_i = freq_df.set_index("gene").loc[gi, "frequency"]
                f_j = freq_df.set_index("gene").loc[gj, "frequency"]
                loser = gi if f_i < f_j else gj
                if loser in keep:
                    keep.remove(loser)
                    dropped.append((loser, gj if loser == gi else gi))
    if dropped:
        print(f"dropped collinear partners (r > 0.999): {dropped}")
    panel = keep

    print(f"\nfinal locked panel ({len(panel)} genes): {panel}")
    if set(panel) != set(config.EXPECTED_PANEL_GENES):
        only_us = set(panel) - set(config.EXPECTED_PANEL_GENES)
        only_them = set(config.EXPECTED_PANEL_GENES) - set(panel)
        print(f"NOTE: panel differs from manuscript reference.")
        print(f"  in this run but not manuscript:  {only_us}")
        print(f"  in manuscript but not this run:  {only_them}")
        print(f"  This is acceptable seed-to-seed variation if the difference is")
        print(f"  among genes whose selection frequency is near the 60% threshold.")

    # Refit on full training set
    y = pheno.set_index("sample_id").loc[X_panel.index, "asd"].astype(int).values
    model = refit_locked_panel(
        X_train=X_panel,
        y_train=y,
        panel_genes=panel,
        lasso_cs=get_lasso_cs(),
        inner_folds=config.INNER_CV_FOLDS,
        lasso_max_iter=config.LASSO_MAX_ITER,
        seed=config.SEED,
    )

    # Apparent (training-fit) AUC and Brier as a sanity check
    Xz = (X_panel - X_panel.mean()) / X_panel.std()
    lp = float(model["intercept"]) + (Xz[panel] * pd.Series(model["coef"], index=panel)).sum(axis=1)
    prob = 1.0 / (1.0 + np.exp(-np.clip(lp.values, -50, 50)))
    model["apparent_auc"] = float(roc_auc_score(y, prob))
    model["apparent_brier"] = float(brier_score_loss(y, prob))
    print(f"\napparent AUC on training: {model['apparent_auc']:.4f}")
    print(f"apparent Brier on training: {model['apparent_brier']:.4f}")

    # Save
    write_json(model, config.LOCKED_MODEL_MAIN)
    print(f"\nwrote {config.LOCKED_MODEL_MAIN}")

    # Also save the reduced 6-gene model (for the GSE42133 Illumina platform)
    reduced_panel = [g for g in panel if g not in config.GPL10558_MISSING_PANEL_GENES]
    if len(reduced_panel) < len(panel):
        X_reduced = X_panel[reduced_panel]
        reduced_model = refit_locked_panel(
            X_train=X_reduced,
            y_train=y,
            panel_genes=reduced_panel,
            lasso_cs=get_lasso_cs(),
            inner_folds=config.INNER_CV_FOLDS,
            lasso_max_iter=config.LASSO_MAX_ITER,
            seed=config.SEED,
        )
        write_json(reduced_model, config.LOCKED_MODEL_REDUCED)
        print(f"wrote {config.LOCKED_MODEL_REDUCED} (panel = {reduced_panel})")

    print("\nstage 05 complete.")


if __name__ == "__main__":
    sys.exit(main() or 0)
