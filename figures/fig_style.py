"""
fig_style.py — Shared style + helpers for iScience submission figures.

Cell Press iScience figure conventions:
  - Sans-serif (Arial / Helvetica)
  - 8 pt body text, 9 pt panel titles, 7 pt tick labels
  - Single column = 85 mm (3.35 in); double column = 114 mm (4.5 in);
    full width = 174 mm (6.85 in)
  - 300 dpi minimum for raster
  - Lines >= 0.5 pt, axis spine 0.75 pt
"""
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, brier_score_loss
from scipy import stats

# ---------- Global rcParams ----------
mpl.rcParams.update({
    "font.family": "DejaVu Sans",   # falls back gracefully; Arial-equivalent metrics
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.75,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "lines.linewidth": 1.2,
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,  # embed TrueType so editors can edit text
    "ps.fonttype": 42,
})

# ---------- Width constants (inches) ----------
WIDTH_SINGLE = 3.35
WIDTH_DOUBLE = 4.5
WIDTH_FULL = 6.85

# ---------- Colour palette (iScience-friendly, colorblind-safe) ----------
COL = {
    "training": "#2c3e50",      # dark slate (training cohort)
    "gpl6244":  "#1f77b4",      # blue (within-series external)
    "gse42133": "#d62728",      # red (independent external)
    "gse25507": "#ff7f0e",      # orange (tissue boundary)
    "asd":      "#c0392b",      # ASD case colour
    "ctrl":     "#27ae60",      # control colour
    "male":     "#3498db",      # male
    "female":   "#e91e63",      # female
    "adapted":  "#1f77b4",
    "strict":   "#9467bd",      # purple (sensitivity)
    "before":   "#7f7f7f",      # pre-recal
    "after":    "#2ca02c",      # post-recal (green)
    "pos":      "#1f77b4",      # positive coefficient
    "neg":      "#d62728",      # negative coefficient
    "diag":     "#bdbdbd",      # diagonal
    "null":     "#bdbdbd",      # permutation null
    "obs":      "#c0392b",      # observed value
}

# ---------- Cell-marker reference (consistent with pipeline.CELL_MARKERS) ----------
CELL_TYPE_ORDER = [
    "Monocytes", "Neutrophils", "Dendritic", "NK cells",
    "CD4 T cells", "CD8 T cells", "B cells", "Erythrocytes",
]

# ---------- Locked model helpers ----------
def load_locked_model(path):
    with open(path) as f:
        m = json.load(f)
    return m


def panel_genes(m):
    return m.get("panel") or m["genes"]


def predict_adapted(model, X):
    """Platform-adapted: z-score with cohort's own mean/SD; returns (P, LP)."""
    panel = panel_genes(model)
    present = [g for g in panel if g in X.columns]
    Xp = X[present].astype(float)
    Z = (Xp - Xp.mean()) / Xp.std()
    if isinstance(model["coefficients"], dict):
        coefs = np.array([model["coefficients"][g] for g in present])
    else:
        coefs = np.array([model["coefficients"][panel.index(g)] for g in present])
    LP = float(model["intercept"]) + Z.values @ coefs
    P = 1.0 / (1.0 + np.exp(-np.clip(LP, -50, 50)))
    return P, LP


def predict_strict(model, X):
    """Strict: z-score with training mean/SD; returns (P, LP)."""
    panel = panel_genes(model)
    present = [g for g in panel if g in X.columns]
    mu = np.array([model["train_mean"][g] for g in present])
    sd = np.array([model["train_sd"][g] for g in present])
    Xp = X[present].astype(float).values
    Z = (Xp - mu) / sd
    if isinstance(model["coefficients"], dict):
        coefs = np.array([model["coefficients"][g] for g in present])
    else:
        coefs = np.array([model["coefficients"][panel.index(g)] for g in present])
    LP = float(model["intercept"]) + Z @ coefs
    P = 1.0 / (1.0 + np.exp(-np.clip(LP, -50, 50)))
    return P, LP


def predict_zerofill_adapted(model, X):
    """7-gene zero-fill with platform-adapted standardisation."""
    panel = panel_genes(model)
    if isinstance(model["coefficients"], dict):
        coefs = np.array([model["coefficients"][g] for g in panel])
    else:
        coefs = np.array(model["coefficients"])
    Z = np.zeros((X.shape[0], len(panel)))
    for j, g in enumerate(panel):
        if g in X.columns:
            v = X[g].astype(float).values
            sd = v.std()
            if sd > 0:
                Z[:, j] = (v - v.mean()) / sd
    LP = float(model["intercept"]) + Z @ coefs
    P = 1.0 / (1.0 + np.exp(-np.clip(LP, -50, 50)))
    return P, LP


def load_cohort(key, work_dir=".."):
    """Returns (expr_samples_x_genes, y, pheno)."""
    expr = pd.read_parquet(f"{work_dir}/expr_gene_{key}.parquet")
    pheno = pd.read_csv(f"{work_dir}/phenotype_{key}.csv")
    # parquet is gene × sample → transpose if needed
    if expr.shape[0] > expr.shape[1] and not str(expr.index[0]).startswith("GSM"):
        expr = expr.T
    sids = pheno["sample_id"].values
    expr = expr.loc[sids]
    y = pheno.set_index("sample_id").loc[expr.index, "label"].astype(int).values
    return expr, y, pheno


# ---------- DeLong CI ----------
def delong_ci(y_true, y_score, alpha=0.05):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n1 = int((y_true == 1).sum())
    n0 = int((y_true == 0).sum())
    if n1 == 0 or n0 == 0:
        return (np.nan, np.nan)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    v10 = np.array([((p > neg).sum() + 0.5 * (p == neg).sum()) / n0 for p in pos])
    v01 = np.array([((pos > n).sum() + 0.5 * (pos == n).sum()) / n1 for n in neg])
    auc = (v10.sum()) / n1
    se = np.sqrt(v10.var(ddof=1) / n1 + v01.var(ddof=1) / n0)
    z = stats.norm.ppf(1 - alpha / 2)
    return max(0.0, auc - z * se), min(1.0, auc + z * se)


# ---------- Calibration curve (binned) ----------
def calibration_bin(y_true, p_pred, n_bins=10, strategy="quantile"):
    y_true = np.asarray(y_true).astype(int)
    p_pred = np.asarray(p_pred, dtype=float)
    if strategy == "quantile":
        edges = np.quantile(p_pred, np.linspace(0, 1, n_bins + 1))
        edges = np.unique(edges)
    else:
        edges = np.linspace(0.0, 1.0, n_bins + 1)
    if len(edges) < 3:
        edges = np.linspace(p_pred.min(), p_pred.max(), n_bins + 1)
    bin_id = np.digitize(p_pred, edges[1:-1])
    pred_means, obs_means, ns = [], [], []
    for b in range(len(edges) - 1):
        mask = bin_id == b
        if mask.sum() == 0:
            continue
        pred_means.append(p_pred[mask].mean())
        obs_means.append(y_true[mask].mean())
        ns.append(mask.sum())
    return np.array(pred_means), np.array(obs_means), np.array(ns)


# ---------- Cell-marker scoring (consistent with pipeline.CELL_MARKERS) ----------
CELL_MARKERS = {
    "Monocytes": ["CD14", "CD163", "FCN1", "VCAN", "S100A8", "S100A9", "LYZ", "CSTA"],
    "Neutrophils": ["FCGR3B", "CSF3R", "CXCR2", "CXCR1", "ELANE", "MPO", "DEFA4", "LCN2"],
    "B cells": ["CD19", "CD79A", "CD79B", "MS4A1", "BLNK", "BANK1", "PAX5", "TCL1A"],
    "CD4 T cells": ["CD4", "IL7R", "TCF7", "LEF1", "CCR7", "FOXP3", "CTLA4"],
    "CD8 T cells": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "NKG7"],
    "NK cells": ["NCAM1", "NKG7", "KLRD1", "KLRC1", "FCGR3A", "GNLY"],
    "Dendritic": ["CLEC10A", "CD1C", "FCER1A", "CLEC4C"],
    "Erythrocytes": ["HBB", "HBA1", "HBA2", "ALAS2", "EPOR"],
}


def cell_marker_score(X, markers):
    """Per-sample mean of z-scored expression across markers present in X."""
    present = [m for m in markers if m in X.columns]
    if not present:
        return None
    Z = (X[present] - X[present].mean()) / X[present].std()
    return Z.mean(axis=1).values


# ---------- Saving helpers ----------
def save_fig(fig, name, out_dir="output"):
    """Save figure as PNG, PDF, TIFF (300 dpi)."""
    from pathlib import Path
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_dir}/{name}.png", dpi=300)
    fig.savefig(f"{out_dir}/{name}.pdf")
    fig.savefig(f"{out_dir}/{name}.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    print(f"  wrote {out_dir}/{name}.{{png,pdf,tiff}}")
