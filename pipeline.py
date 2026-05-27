"""Shared analytic functions used across the pipeline.

Conventions:
- expression matrices are stored as `pandas.DataFrame` with samples as rows
  and genes (or probes) as columns;
- phenotype tables are stored as `pandas.DataFrame` with a `sample_id` column,
  an `asd` column (1 = ASD, 0 = typically developing), and optional `sex` /
  `age` columns;
- all log2 transformations are upstream of this module — these functions
  assume their inputs are already on the log2 scale.
"""
from __future__ import annotations

import gzip
import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, brier_score_loss
from joblib import Parallel, delayed

import config


# ============================================================================
# I/O helpers
# ============================================================================

def write_json(obj, path: Path) -> None:
    """Write a JSON file with deterministic key ordering and two-space indent."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=_json_default)


def read_json(path: Path):
    with open(path) as f:
        return json.load(f)


def _json_default(o):
    """JSON default serialiser for numpy scalars and numpy arrays."""
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"Object of type {type(o)} is not JSON serialisable")


# ============================================================================
# Preprocessing: probe → gene IQR-max aggregation
# ============================================================================

def iqr_max_aggregate(expr_probe: pd.DataFrame,
                      probe_to_gene: pd.Series) -> pd.DataFrame:
    """Collapse a probe-level expression matrix to gene level.

    For each gene with multiple mapped probes, retain the probe with the
    largest interquartile range (IQR) across samples. Drops probes whose
    gene symbol is missing or blank.

    Parameters
    ----------
    expr_probe : DataFrame, samples (rows) × probes (cols)
    probe_to_gene : Series indexed by probe id, values are gene symbols

    Returns
    -------
    DataFrame, samples (rows) × genes (cols), sorted by gene symbol
    """
    # Compute IQR per probe (over samples → axis=0)
    q75 = expr_probe.quantile(0.75, axis=0)
    q25 = expr_probe.quantile(0.25, axis=0)
    iqr = (q75 - q25)

    # Map probe → gene; drop probes with no gene assignment
    genes = expr_probe.columns.map(probe_to_gene)
    keep = ~(pd.isna(genes) | (genes.astype(str).str.strip() == ""))

    aux = pd.DataFrame({
        "probe": expr_probe.columns[keep],
        "gene": genes[keep],
        "iqr": iqr.values[keep],
    })

    # For each gene, keep the probe with highest IQR (ties → first)
    aux = aux.sort_values(["gene", "iqr"], ascending=[True, False])
    aux = aux.drop_duplicates("gene", keep="first")

    keep_probes = aux["probe"].tolist()
    gene_for_probe = dict(zip(aux["probe"], aux["gene"]))
    out = expr_probe[keep_probes].copy()
    out.columns = [gene_for_probe[p] for p in keep_probes]
    out = out.sort_index(axis=1)
    return out


# ============================================================================
# Standardisation
# ============================================================================

def zscore_using(stats_mean: pd.Series, stats_sd: pd.Series,
                 X: pd.DataFrame) -> pd.DataFrame:
    """Return (X − μ) / σ using the supplied per-gene mean and SD.

    Genes in `X` that are missing from `stats_mean` are dropped from the
    output; the caller is responsible for handling missing panel genes
    (e.g. the GPL10558 missing-ERVK3-2 case in 09_external_GSE42133.py).
    """
    common = [g for g in stats_mean.index if g in X.columns]
    Z = (X[common] - stats_mean[common]) / stats_sd[common]
    return Z


def cohort_internal_zscore(X: pd.DataFrame) -> pd.DataFrame:
    """Per-column z-score using cohort-internal mean and SD."""
    return (X - X.mean()) / X.std()


# ============================================================================
# Linear predictor and probability
# ============================================================================

def linear_predictor(Z: pd.DataFrame, coefs: pd.Series, intercept: float) -> pd.Series:
    """LP = β₀ + Σ βᵢ · zᵢ.

    Only the intersection of `Z.columns` and `coefs.index` is used; missing
    panel genes contribute 0 to the linear predictor (caller controls
    whether this matches the intended sensitivity analysis).
    """
    common = [g for g in coefs.index if g in Z.columns]
    return float(intercept) + (Z[common] * coefs[common]).sum(axis=1)


def sigmoid(x):
    x = np.clip(x, -50, 50)  # numerical stability
    return 1.0 / (1.0 + np.exp(-x))


def predict_proba_locked(model: dict, X: pd.DataFrame,
                          mode: str = "platform_adapted") -> pd.Series:
    """Apply a locked model to an external cohort.

    Parameters
    ----------
    model : dict with keys 'genes', 'coef', 'intercept', 'train_mean',
            'train_sd' (as produced by lock_and_refit_panel)
    X : DataFrame, samples × genes (log2 scale, gene-aggregated)
    mode : one of:
        'platform_adapted'  — z-score with this cohort's own mean / SD
        'strict_train'      — z-score with the training-set mean / SD
        'zero_fill'         — strict_train but z = 0 for genes missing in X
    """
    coefs = pd.Series(model["coef"], index=model["genes"])
    intercept = float(model["intercept"])

    if mode == "platform_adapted":
        # cohort-internal standardisation across this cohort
        in_X = [g for g in model["genes"] if g in X.columns]
        if not in_X:
            raise ValueError("None of the panel genes are present in X.")
        X_sub = X[in_X]
        Z = (X_sub - X_sub.mean()) / X_sub.std()
    elif mode == "strict_train":
        mu = pd.Series(model["train_mean"], index=model["genes"])
        sd = pd.Series(model["train_sd"], index=model["genes"])
        in_X = [g for g in model["genes"] if g in X.columns]
        if not in_X:
            raise ValueError("None of the panel genes are present in X.")
        Z = (X[in_X] - mu[in_X]) / sd[in_X]
    elif mode == "zero_fill":
        mu = pd.Series(model["train_mean"], index=model["genes"])
        sd = pd.Series(model["train_sd"], index=model["genes"])
        Z = pd.DataFrame(0.0, index=X.index, columns=model["genes"])
        present = [g for g in model["genes"] if g in X.columns]
        Z.loc[:, present] = ((X[present] - mu[present]) / sd[present]).values
    else:
        raise ValueError(f"unknown mode {mode!r}")

    lp = linear_predictor(Z, coefs, intercept)
    p = sigmoid(lp.values)
    return pd.Series(p, index=X.index, name="prob")


# ============================================================================
# Statistical helpers
# ============================================================================

def auc(y_true, y_score) -> float:
    return float(roc_auc_score(y_true, y_score))


def brier(y_true, y_prob) -> float:
    return float(brier_score_loss(y_true, y_prob))


def delong_ci(y_true: np.ndarray, y_score: np.ndarray, alpha: float = 0.05):
    """95% CI for the AUC using DeLong's method.

    Implementation follows Sun & Xu (2014) reformulation; assumes binary
    outcome and one classifier. Returns (lower, upper).
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n1 = int((y_true == 1).sum())
    n0 = int((y_true == 0).sum())
    if n1 == 0 or n0 == 0:
        return (float("nan"), float("nan"))

    # Mid-rank for ties
    order = np.argsort(y_score, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    # Handle ties using fractional ranks (average rank for tied values)
    s = y_score[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1

    pos_ranks = ranks[y_true == 1]
    neg_ranks = ranks[y_true == 0]
    auc_val = (pos_ranks.sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)

    # Structural components (Sun-Xu)
    # V10[i] = (1 / n0) * sum_j I(yp_i > yn_j) + 0.5 * I(yp_i == yn_j)
    pos_scores = y_score[y_true == 1]
    neg_scores = y_score[y_true == 0]

    v10 = np.empty(n1)
    v01 = np.empty(n0)
    for k in range(n1):
        v10[k] = ((pos_scores[k] > neg_scores).sum()
                  + 0.5 * (pos_scores[k] == neg_scores).sum()) / n0
    for k in range(n0):
        v01[k] = ((pos_scores > neg_scores[k]).sum()
                  + 0.5 * (pos_scores == neg_scores[k]).sum()) / n1

    s10 = v10.var(ddof=1) / n1
    s01 = v01.var(ddof=1) / n0
    se = math.sqrt(s10 + s01)
    z = stats.norm.ppf(1 - alpha / 2)
    return (max(0.0, auc_val - z * se), min(1.0, auc_val + z * se))


def permutation_pvalue(y_true: np.ndarray, y_score: np.ndarray,
                       n_perms: int = 10000, seed: int = 42,
                       two_sided: bool = True) -> float:
    """Empirical p-value: P(perm AUC ≥ observed) under shuffled labels.

    For two-sided, returns 2 × min(upper-tail, lower-tail), capped at 1.
    """
    rng = np.random.RandomState(seed)
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    observed = roc_auc_score(y_true, y_score)

    n_upper = 0
    n_lower = 0
    y_perm = y_true.copy()
    for _ in range(n_perms):
        rng.shuffle(y_perm)
        a = roc_auc_score(y_perm, y_score)
        if a >= observed:
            n_upper += 1
        if a <= observed:
            n_lower += 1

    upper = (n_upper + 1) / (n_perms + 1)
    lower = (n_lower + 1) / (n_perms + 1)
    if two_sided:
        return min(1.0, 2.0 * min(upper, lower))
    return upper


def logistic_recalibration(y_true: np.ndarray, lp: np.ndarray):
    """Fit p = sigmoid(α + β · LP) on the supplied data.

    Returns (alpha, beta).
    """
    # Build a 1-feature logistic regression (no penalty) on the linear predictor
    lp = np.asarray(lp, dtype=float).reshape(-1, 1)
    y_true = np.asarray(y_true).astype(int)
    clf = LogisticRegression(C=1e8, solver="lbfgs", max_iter=1000)
    clf.fit(lp, y_true)
    beta = float(clf.coef_[0, 0])
    alpha = float(clf.intercept_[0])
    return alpha, beta


# ============================================================================
# Nested cross-validation with t-test pre-screen and LASSO
# ============================================================================

def _ttest_top_k(X_train: pd.DataFrame, y_train: np.ndarray,
                 top_k: int) -> list[str]:
    """Welch's two-sample t-test, return top-k gene names by ascending p-value."""
    grp1 = X_train.values[y_train == 1]
    grp0 = X_train.values[y_train == 0]
    t, p = stats.ttest_ind(grp1, grp0, equal_var=False, axis=0)
    p = np.nan_to_num(p, nan=1.0)
    top_idx = np.argsort(p)[:top_k]
    return [X_train.columns[i] for i in top_idx]


def _outer_fold_worker(args):
    """Single outer-fold worker for joblib parallelism.

    Returns dict with: repeat, fold, auc, selected_genes (list of str)
    """
    (rep, fold_idx, train_idx, test_idx, X_values, y_values,
     gene_names, top_k, lasso_cs, lasso_max_iter, inner_folds, seed) = args

    X_train = X_values[train_idx]
    X_test = X_values[test_idx]
    y_train = y_values[train_idx]
    y_test = y_values[test_idx]

    # t-test pre-screen on training subset only
    t, p = stats.ttest_ind(X_train[y_train == 1],
                           X_train[y_train == 0],
                           equal_var=False, axis=0)
    p = np.nan_to_num(p, nan=1.0)
    top_idx = np.argsort(p)[:top_k]
    Xtr = X_train[:, top_idx]
    Xte = X_test[:, top_idx]

    # Standardise using training subset only
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    Xtr_z = (Xtr - mu) / sd
    Xte_z = (Xte - mu) / sd

    # Inner CV LASSO
    inner_cv = StratifiedKFold(n_splits=inner_folds, shuffle=True,
                               random_state=int(seed))
    clf = LogisticRegressionCV(
        Cs=lasso_cs, penalty="l1", solver="liblinear",
        cv=inner_cv, scoring="roc_auc",
        max_iter=lasso_max_iter, random_state=int(seed),
    )
    clf.fit(Xtr_z, y_train)

    prob = clf.predict_proba(Xte_z)[:, 1]
    fold_auc = float(roc_auc_score(y_test, prob))

    coef = clf.coef_[0]
    selected = [gene_names[top_idx[j]] for j in range(len(coef)) if coef[j] != 0.0]
    return {"repeat": int(rep), "fold": int(fold_idx), "auc": fold_auc,
            "selected_genes": selected,
            "n_selected": int(len(selected))}


def nested_cv_lasso_stability(X: pd.DataFrame, y: pd.Series,
                              n_outer: int, n_repeats: int,
                              top_k: int, inner_folds: int,
                              lasso_cs: np.ndarray, lasso_max_iter: int,
                              seed: int, n_jobs: int = -1) -> list[dict]:
    """Run the full nested-CV pipeline; return list of per-fold result dicts."""
    rng = np.random.RandomState(seed)
    X_values = X.values
    y_values = np.asarray(y).astype(int)
    gene_names = list(X.columns)

    # Build the master task list (one tuple per outer fold across all repeats)
    tasks = []
    for rep in range(n_repeats):
        outer_seed = int(rng.randint(2**31 - 1))
        cv = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=outer_seed)
        for fold_idx, (tr, te) in enumerate(cv.split(X_values, y_values)):
            tasks.append((
                rep, fold_idx, tr, te, X_values, y_values, gene_names,
                top_k, lasso_cs, lasso_max_iter, inner_folds, outer_seed + fold_idx,
            ))

    print(f"  nested CV: {len(tasks)} outer folds, n_jobs={n_jobs}", flush=True)
    results = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_outer_fold_worker)(t) for t in tasks
    )
    return results


def lock_panel_from_nested_results(results: list[dict],
                                    threshold: float) -> tuple[list[str], dict]:
    """Return (panel, freq_map) given nested-CV per-fold results."""
    n_folds = len(results)
    counts = Counter()
    for r in results:
        for g in r["selected_genes"]:
            counts[g] += 1
    freq = {g: c / n_folds for g, c in counts.items()}
    panel = sorted([g for g, f in freq.items() if f >= threshold],
                   key=lambda g: -freq[g])
    return panel, freq


def refit_locked_panel(X_train: pd.DataFrame, y_train: np.ndarray,
                       panel_genes: list[str],
                       lasso_cs: np.ndarray, inner_folds: int,
                       lasso_max_iter: int, seed: int) -> dict:
    """Refit LASSO on the locked panel only, on the full training set.

    Returns dict with keys: genes, coef, intercept, train_mean, train_sd.
    """
    missing = [g for g in panel_genes if g not in X_train.columns]
    if missing:
        raise ValueError(f"panel genes missing from training expr: {missing}")

    Xp = X_train[panel_genes].copy()
    mu = Xp.mean()
    sd = Xp.std()
    sd[sd == 0] = 1.0
    Xz = (Xp - mu) / sd

    inner_cv = StratifiedKFold(n_splits=inner_folds, shuffle=True, random_state=seed)
    clf = LogisticRegressionCV(
        Cs=lasso_cs, penalty="l1", solver="liblinear",
        cv=inner_cv, scoring="roc_auc",
        max_iter=lasso_max_iter, random_state=seed,
    )
    clf.fit(Xz.values, y_train)

    coef = clf.coef_[0]
    intercept = float(clf.intercept_[0])
    return {
        "genes": panel_genes,
        "coef": [float(c) for c in coef],
        "intercept": intercept,
        "train_mean": {g: float(mu[g]) for g in panel_genes},
        "train_sd": {g: float(sd[g]) for g in panel_genes},
        "C_chosen": float(clf.C_[0]),
    }


# ============================================================================
# 0.632 bootstrap optimism correction (panel held fixed)
# ============================================================================

def bootstrap_optimism(X_train: pd.DataFrame, y_train: np.ndarray,
                       panel_genes: list[str],
                       lasso_cs: np.ndarray, inner_folds: int,
                       lasso_max_iter: int, B: int, seed: int) -> dict:
    """0.632 optimism correction: apparent AUC − mean(optimism).

    Holds the gene panel fixed; only refits coefficients per bootstrap sample.
    """
    rng = np.random.RandomState(seed)
    n = len(y_train)
    Xp = X_train[panel_genes].values
    y = np.asarray(y_train).astype(int)

    # Apparent AUC: refit on full data, evaluate on full data
    mu = Xp.mean(axis=0)
    sd = Xp.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    Xz = (Xp - mu) / sd
    inner_cv = StratifiedKFold(n_splits=inner_folds, shuffle=True, random_state=seed)
    clf = LogisticRegressionCV(Cs=lasso_cs, penalty="l1", solver="liblinear",
                                cv=inner_cv, scoring="roc_auc",
                                max_iter=lasso_max_iter, random_state=seed)
    clf.fit(Xz, y)
    apparent = float(roc_auc_score(y, clf.predict_proba(Xz)[:, 1]))

    optimisms = []
    for b in range(B):
        idx = rng.choice(n, size=n, replace=True)
        X_b = Xp[idx]
        y_b = y[idx]
        # Bootstrap standardisation
        mu_b = X_b.mean(axis=0)
        sd_b = X_b.std(axis=0, ddof=1)
        sd_b[sd_b == 0] = 1.0
        X_b_z = (X_b - mu_b) / sd_b
        # Apply same standardisation to original data
        X_orig_z = (Xp - mu_b) / sd_b

        cv_b = StratifiedKFold(n_splits=inner_folds, shuffle=True,
                               random_state=seed + b + 1)
        try:
            clf_b = LogisticRegressionCV(Cs=lasso_cs, penalty="l1", solver="liblinear",
                                          cv=cv_b, scoring="roc_auc",
                                          max_iter=lasso_max_iter,
                                          random_state=seed + b + 1)
            clf_b.fit(X_b_z, y_b)
            auc_train = float(roc_auc_score(y_b, clf_b.predict_proba(X_b_z)[:, 1]))
            auc_test = float(roc_auc_score(y, clf_b.predict_proba(X_orig_z)[:, 1]))
            optimisms.append(auc_train - auc_test)
        except ValueError:
            # Singular bootstrap sample (all one class) — skip
            continue

        if (b + 1) % 100 == 0:
            print(f"  bootstrap {b + 1}/{B} done", flush=True)

    mean_optimism = float(np.mean(optimisms))
    corrected = apparent - mean_optimism
    return {
        "apparent_auc": apparent,
        "mean_optimism": mean_optimism,
        "optimism_corrected_auc": corrected,
        "B_used": len(optimisms),
        "B_requested": B,
    }


# ============================================================================
# Cell-type marker reference panels
# ============================================================================
# Compiled from CIBERSORT LM22 (Newman et al. 2015) and xCell (Aran et al. 2017)
# Each list is the consensus marker set used in this study.

CELL_MARKERS = {
    "monocytes": ["CD14", "CD163", "FCN1", "VCAN", "S100A8", "S100A9", "LYZ", "CSTA"],
    "neutrophils": ["FCGR3B", "CSF3R", "CXCR2", "CXCR1", "ELANE", "MPO", "DEFA4", "LCN2"],
    "B_cells": ["CD19", "CD79A", "CD79B", "MS4A1", "BLNK", "BANK1", "PAX5", "TCL1A"],
    "CD4_T_cells": ["CD4", "IL7R", "TCF7", "LEF1", "CCR7", "FOXP3", "CTLA4"],
    "CD8_T_cells": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "NKG7"],
    "NK_cells": ["NCAM1", "NKG7", "KLRD1", "KLRC1", "FCGR3A", "GNLY"],
    "dendritic_cells": ["CLEC10A", "CD1C", "FCER1A", "CLEC4C"],
    "erythrocytes": ["HBB", "HBA1", "HBA2", "ALAS2", "EPOR"],
}


def cell_type_marker_score(X: pd.DataFrame, markers: list[str]) -> pd.Series:
    """Per-sample mean of z-scored expression across the supplied markers.

    Markers absent from `X.columns` are silently skipped.
    """
    present = [m for m in markers if m in X.columns]
    if not present:
        return pd.Series(np.nan, index=X.index)
    Z = (X[present] - X[present].mean()) / X[present].std()
    return Z.mean(axis=1)


def cell_marker_correlations(X: pd.DataFrame, lp7: pd.Series) -> dict:
    """For each cell type, return Pearson r between LP7 and the marker score."""
    out = {}
    for ct, markers in CELL_MARKERS.items():
        score = cell_type_marker_score(X, markers)
        valid = lp7.index.intersection(score.index)
        if len(valid) < 3 or score.loc[valid].isna().all():
            out[ct] = {"r": float("nan"), "p": float("nan"), "n": 0,
                       "n_markers": 0}
            continue
        r, p = stats.pearsonr(lp7.loc[valid].values, score.loc[valid].values)
        out[ct] = {"r": float(r), "p": float(p), "n": int(len(valid)),
                   "n_markers": int(len([m for m in markers if m in X.columns]))}
    return out


# ============================================================================
# Lasso grid (used everywhere; defined once for consistency)
# ============================================================================

def get_lasso_cs() -> np.ndarray:
    return np.logspace(config.LASSO_LOG10_C_MIN,
                       config.LASSO_LOG10_C_MAX,
                       config.LASSO_N_CS)
