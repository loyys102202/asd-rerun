"""
make_sfig2.py — Supplementary Figure 2: LASSO regularisation path (v2).

Improvements over v1:
- 10× repeated 5-fold CV for smoother CV-AUC curve
- Locked-panel coefficient paths drawn at high zorder over background
- λ_min and stability-selection-active region marked
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

import fig_style as fs

# ---------- Load ----------
model = fs.load_locked_model("../locked_model_main_7gene.json")
panel = model["panel"]
expr, y, _ = fs.load_cohort("GPL570")

# Standardise
mu = expr.mean(); sd = expr.std()
Z = ((expr - mu) / sd).astype(float).replace([np.inf, -np.inf], 0.0).fillna(0.0)
X_all = Z.values
y = y.astype(int)
gene_names = np.array(expr.columns)

# ---------- λ grid (a bit tighter than v1) ----------
log_lambdas = np.linspace(-4.5, 0.5, 30)
Cs = 1.0 / (10.0 ** log_lambdas)

# ---------- Coefficient paths (fit on full training set at each λ) ----------
coef_path = np.zeros((X_all.shape[1], len(Cs)))
for i, C in enumerate(Cs):
    clf = LogisticRegression(penalty="l1", solver="liblinear",
                              C=C, max_iter=5000, random_state=42)
    clf.fit(X_all, y)
    coef_path[:, i] = clf.coef_[0]

ever_entered = np.where(np.abs(coef_path).max(axis=1) > 1e-6)[0]
panel_idx_set = set([list(gene_names).index(g) for g in panel if g in gene_names])
others_idx = [i for i in ever_entered if i not in panel_idx_set]
print(f"  Total genes ever entering: {len(ever_entered)}")

# ---------- Repeated 5-fold CV AUC per λ (10 repeats × 5 folds = 50 evaluations) ----------
N_REPEATS = 5
auc_means = np.zeros(len(Cs))
auc_sds = np.zeros(len(Cs))
for i, C in enumerate(Cs):
    aucs = []
    for rep in range(N_REPEATS):
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=rep * 31 + 7)
        for tr_idx, va_idx in cv.split(X_all, y):
            clf = LogisticRegression(penalty="l1", solver="liblinear",
                                      C=C, max_iter=5000, random_state=42)
            clf.fit(X_all[tr_idx], y[tr_idx])
            P = clf.predict_proba(X_all[va_idx])[:, 1]
            aucs.append(roc_auc_score(y[va_idx], P))
    auc_means[i] = np.mean(aucs)
    auc_sds[i] = np.std(aucs)
    if i % 10 == 0:
        print(f"    λ={10**log_lambdas[i]:.3g}  AUC mean={auc_means[i]:.3f}  sd={auc_sds[i]:.3f}")

idx_min = int(np.argmax(auc_means))
log_lambda_min = log_lambdas[idx_min]
auc_min = auc_means[idx_min]
sd_at_min = auc_sds[idx_min]
target = auc_min - sd_at_min
right = np.where(auc_means[idx_min:] >= target)[0]
idx_1se = idx_min + (right[-1] if len(right) > 0 else 0)
log_lambda_1se = log_lambdas[idx_1se]

print(f"  λ_min = 10^{log_lambda_min:.2f}, CV-AUC = {auc_min:.3f} ± {sd_at_min:.3f}")
print(f"  λ_1se = 10^{log_lambda_1se:.2f}, CV-AUC = {auc_means[idx_1se]:.3f}")

# ---------- Figure ----------
fig = plt.figure(figsize=(fs.WIDTH_FULL, 4.0))
gs = gridspec.GridSpec(nrows=1, ncols=2,
                       left=0.07, right=0.97, top=0.91, bottom=0.27,
                       wspace=0.40)

panel_colors = {
    "KYNU": "#1f77b4", "CMYA5": "#2ca02c", "MAPK8IP1": "#9467bd",
    "CES1": "#17becf", "KIAA0087": "#d62728", "POU2AF1": "#ff7f0e",
    "ERVK3-2": "#8c564b",
}

# ===== Panel A =====
axA = fig.add_subplot(gs[0, 0])
# Background — other genes at low zorder
for i in others_idx:
    axA.plot(log_lambdas, coef_path[i, :], color="#cccccc",
              linewidth=0.4, alpha=0.5, zorder=1)
# Locked-panel genes — bold lines at high zorder
for g in panel:
    if g in gene_names:
        gi = list(gene_names).index(g)
        axA.plot(log_lambdas, coef_path[gi, :], color=panel_colors[g],
                  linewidth=1.8, label=g, zorder=5)

axA.axvline(log_lambda_min, color="black", linewidth=1.0, linestyle="--", zorder=4)
axA.axhline(0, color="black", linewidth=0.4, zorder=2)
axA.set_xlabel("log$_{10}$(λ)", fontsize=8.5)
axA.set_ylabel("Logit-scale coefficient", fontsize=8.5)
axA.set_title("A. LASSO coefficient paths on the training set",
              loc="left", fontweight="bold", fontsize=9)
axA.set_xlim(log_lambdas.min(), log_lambdas.max())

legend_handles = [mlines.Line2D([0],[0], color=panel_colors[g], linewidth=1.8, label=g)
                  for g in panel]
legend_handles.append(mlines.Line2D([0],[0], color="#cccccc", linewidth=0.4,
                                      label="other genes entering at some λ"))
legend_handles.append(mlines.Line2D([0],[0], color="black", linewidth=1.0,
                                      linestyle="--",
                                      label=f"λ$_{{min}}$ = 10$^{{{log_lambda_min:.2f}}}$"))
axA.legend(handles=legend_handles, loc="upper center",
            bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False, fontsize=6.3,
            handlelength=1.8, handletextpad=0.5, columnspacing=0.9)

# ===== Panel B =====
axB = fig.add_subplot(gs[0, 1])
axB.fill_between(log_lambdas, auc_means - auc_sds, auc_means + auc_sds,
                  color="#1f77b4", alpha=0.20, zorder=1)
axB.plot(log_lambdas, auc_means, color="#1f77b4", linewidth=1.6, zorder=3)
axB.axvline(log_lambda_min, color="black", linewidth=1.0, linestyle="--", zorder=4)
axB.axvline(log_lambda_1se, color="#7f7f7f", linewidth=1.0, linestyle=":", zorder=4)

axB.set_xlabel("log$_{10}$(λ)", fontsize=8.5)
axB.set_ylabel("Repeated 5-fold CV AUC", fontsize=8.5)
axB.set_title("B. CV-AUC as a function of regularisation",
              loc="left", fontweight="bold", fontsize=9)
axB.set_xlim(log_lambdas.min(), log_lambdas.max())
axB.set_ylim(0.45, 1.0)

legend_b = [
    mlines.Line2D([0],[0], color="#1f77b4", linewidth=1.6,
                  label=f"5-fold CV AUC  ({N_REPEATS} repeats × 5 folds)"),
    mlines.Line2D([0],[0], color="#1f77b4", linewidth=8, alpha=0.20,
                  label="± 1 SD across all folds"),
    mlines.Line2D([0],[0], color="black", linewidth=1.0, linestyle="--",
                  label=f"λ$_{{min}}$ (10$^{{{log_lambda_min:.2f}}}$): AUC = {auc_min:.3f}"),
    mlines.Line2D([0],[0], color="#7f7f7f", linewidth=1.0, linestyle=":",
                  label=f"λ$_{{1se}}$ (10$^{{{log_lambda_1se:.2f}}}$): AUC = {auc_means[idx_1se]:.3f}"),
]
axB.legend(handles=legend_b, loc="upper center",
            bbox_to_anchor=(0.5, -0.18), ncol=1, frameon=False, fontsize=6.5,
            handlelength=2.0)

# Save
out_dir = "output"
os.makedirs(out_dir, exist_ok=True)
fig.savefig(f"{out_dir}/Supplementary_Figure_2.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{out_dir}/Supplementary_Figure_2.pdf", bbox_inches="tight")
fig.savefig(f"{out_dir}/Supplementary_Figure_2.tiff", dpi=300,
             pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")
plt.close(fig)
print("Supplementary Figure 2 v2 saved (PNG + PDF + TIFF).")
