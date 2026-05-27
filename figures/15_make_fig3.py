"""
make_fig3.py v6 — Figure 3: internal performance on GSE18123-GPL570.

Was Figure 2 in the previous 5-figure layout.  Increased inter-column
spacing so each panel's caption/legend strip stays well within its column.
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from sklearn.metrics import roc_curve, roc_auc_score, brier_score_loss

import fig_style as fs

model = fs.load_locked_model("../locked_model_main_7gene.json")
expr, y, pheno = fs.load_cohort("GPL570")
P, LP = fs.predict_strict(model, expr)

auc_by_rep = pd.read_csv("../nestedcv_auc_by_rep.csv")
rep_aucs = auc_by_rep["auc_rep_mean"].values

with open("../robustness_checks.json") as f:
    rc = json.load(f)
with open("../calibration.json") as f:
    cal = json.load(f)

apparent_auc = float(model["apparent_auc"])
nested_mean = float(rep_aucs.mean())
nested_lo = float(np.percentile(rep_aucs, 2.5))
nested_hi = float(np.percentile(rep_aucs, 97.5))
boot_corr = float(rc["bootstrap_optimism_corrected_auc"])
optimism = float(rc["bootstrap_optimism"])
perm_null_mean = float(rc["perm_train_null_mean"])
perm_null_95 = rc["perm_train_null_95"]

fig = plt.figure(figsize=(fs.WIDTH_FULL, 3.6))
gs = gridspec.GridSpec(
    nrows=1, ncols=3,
    left=0.07, right=0.98, top=0.88, bottom=0.30,
    wspace=0.60,
)

# ===== Panel A =====
axA = fig.add_subplot(gs[0, 0])
fpr, tpr, _ = roc_curve(y, P)
axA.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axA.plot(fpr, tpr, color=fs.COL["training"], linewidth=1.6)
axA.set_xlabel("False positive rate")
axA.set_ylabel("True positive rate")
axA.set_title("A. Training-set ROC", loc="left", fontweight="bold", fontsize=9)
axA.set_xlim(-0.02, 1.02); axA.set_ylim(-0.02, 1.02)
axA.set_xticks([0, 0.5, 1.0]); axA.set_yticks([0, 0.5, 1.0])
axA.text(0.5, -0.32,
         f"Apparent AUC  =  {apparent_auc:.3f}\n"
         f"Nested CV    =  {nested_mean:.3f}  [{nested_lo:.3f}, {nested_hi:.3f}]\n"
         f"Bootstrap .632=  {boot_corr:.3f}  (optimism {optimism:.3f})\n"
         f"Permutation p <  0.001",
         transform=axA.transAxes, ha="center", va="top",
         fontsize=6.5, family="monospace",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f5f5f5",
                   edgecolor="#cccccc", linewidth=0.5))

# ===== Panel B =====
axB = fig.add_subplot(gs[0, 1])
axB.hist(rep_aucs, bins=18, color=fs.COL["training"], alpha=0.55,
         edgecolor="black", linewidth=0.4)
axB.axvspan(perm_null_95[0], perm_null_95[1], color=fs.COL["null"], alpha=0.22)
axB.axvline(apparent_auc, color=fs.COL["obs"], linestyle="--", linewidth=1.0)
axB.axvline(boot_corr, color=fs.COL["after"], linestyle="-", linewidth=1.0)
axB.axvline(nested_mean, color=fs.COL["training"], linestyle="-", linewidth=1.4)
axB.set_xlabel("AUC")
axB.set_ylabel("Frequency (per-repeat means)")
axB.set_title("B. AUC distribution", loc="left", fontweight="bold", fontsize=9)
axB.set_xlim(0.30, 1.00)

# Short labels — keep legend narrow enough to stay in this column
hist_p = mpatches.Patch(facecolor=fs.COL["training"], alpha=0.55,
                         edgecolor="black", label="100 nested-CV repeats")
null_p = mpatches.Patch(facecolor=fs.COL["null"], alpha=0.22,
                         label=f"Perm. null 95% (mean {perm_null_mean:.2f})")
app_l = mlines.Line2D([0],[0], color=fs.COL["obs"], linestyle="--", linewidth=1.0,
                       label=f"Apparent ({apparent_auc:.2f})")
boo_l = mlines.Line2D([0],[0], color=fs.COL["after"], linewidth=1.0,
                       label=f"Bootstrap ({boot_corr:.2f})")
nc_l  = mlines.Line2D([0],[0], color=fs.COL["training"], linewidth=1.4,
                       label=f"Nested CV ({nested_mean:.2f})")
axB.legend(handles=[hist_p, null_p, app_l, boo_l, nc_l],
           loc="upper center", bbox_to_anchor=(0.5, -0.24),
           ncol=1, frameon=False, fontsize=6.3,
           handlelength=1.8, handletextpad=0.5)

# ===== Panel C =====
axC = fig.add_subplot(gs[0, 2])
pred_means, obs_means, ns = fs.calibration_bin(y, P, n_bins=10, strategy="quantile")
axC.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axC.plot(pred_means, obs_means, "o-", color=fs.COL["training"], linewidth=1.4,
         markersize=4.5)
axC.set_xlabel("Predicted probability of ASD")
axC.set_ylabel("Observed proportion ASD")
axC.set_title("C. Training calibration", loc="left", fontweight="bold", fontsize=9)
axC.set_xlim(-0.02, 1.02); axC.set_ylim(-0.02, 1.02)
axC.set_xticks([0, 0.5, 1.0]); axC.set_yticks([0, 0.5, 1.0])

brier_val = brier_score_loss(y, P)
cal_int = cal["training_apparent"]["cal_intercept"]
cal_slope = cal["training_apparent"]["cal_slope"]
axC.text(0.5, -0.32,
         f"Brier  =  {brier_val:.3f}\n"
         f"α = {cal_int:+.2f},  β = {cal_slope:+.2f}\n"
         f"Dashed gray: ideal",
         transform=axC.transAxes, ha="center", va="top", fontsize=6.5,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f5f5f5",
                   edgecolor="#cccccc", linewidth=0.5))

fs.save_fig(fig, "Figure_3")
plt.close(fig)
print("Figure 3 done.")
