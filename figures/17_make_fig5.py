"""make_fig5.py v6 — Figure 5: GSE42133 + GSE25507 ROCs (was Fig 4 in v5)."""
import json, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
import fig_style as fs

model7 = fs.load_locked_model("../locked_model_main_7gene.json")
model6 = fs.load_locked_model("../locked_model_reduced_6gene.json")
expr_42133, y_42133, _ = fs.load_cohort("GSE42133")
expr_25507, y_25507, _ = fs.load_cohort("GSE25507")
P_42133, _ = fs.predict_adapted(model6, expr_42133)
P_25507, _ = fs.predict_adapted(model7, expr_25507)
with open("../GSE25507_exploration.json") as f:
    g25 = json.load(f)

fig = plt.figure(figsize=(fs.WIDTH_DOUBLE * 1.5, 3.8))
gs = gridspec.GridSpec(nrows=1, ncols=2, left=0.09, right=0.97,
                       top=0.92, bottom=0.32, wspace=0.40)

# Panel A
axA = fig.add_subplot(gs[0, 0])
fpr, tpr, _ = roc_curve(y_42133, P_42133)
auc = roc_auc_score(y_42133, P_42133); ci = fs.delong_ci(y_42133, P_42133)
axA.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axA.plot(fpr, tpr, color=fs.COL["gse42133"], linewidth=1.5)
axA.set_xlabel("False positive rate"); axA.set_ylabel("True positive rate")
axA.set_title("A. GSE42133 (independent external)",
              loc="left", fontweight="bold", fontsize=9)
axA.set_xlim(-0.02, 1.02); axA.set_ylim(-0.02, 1.02)
axA.set_xticks([0, 0.5, 1.0]); axA.set_yticks([0, 0.5, 1.0])
h_42 = mlines.Line2D([0],[0], color=fs.COL["gse42133"], linewidth=1.5,
                       label=f"GSE42133 (n=147, all male, Illumina HT-12 v4):\n"
                             f"AUC = {auc:.3f}  [{ci[0]:.3f}, {ci[1]:.3f}],  "
                             f"permutation p = 0.12")
h_diag = mlines.Line2D([0],[0], color=fs.COL["diag"], linestyle="--",
                        linewidth=0.7, label="Chance")
axA.legend(handles=[h_42, h_diag], loc="upper center",
            bbox_to_anchor=(0.5, -0.25), ncol=1, frameon=False,
            fontsize=6.8, handlelength=2.0)

# Panel B
axB = fig.add_subplot(gs[0, 1])
fpr, tpr, _ = roc_curve(y_25507, P_25507)
auc = roc_auc_score(y_25507, P_25507); ci = fs.delong_ci(y_25507, P_25507)
axB.plot(fpr, tpr, color=fs.COL["gse25507"], linewidth=1.5)
panel = model7["panel"]
present = [g for g in panel if g in expr_25507.columns]
Xp = expr_25507[present].astype(float); Z = (Xp - Xp.mean()) / Xp.std()
inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
clf = LogisticRegressionCV(Cs=np.logspace(-3, 1, 30), penalty="l1",
                            solver="liblinear", cv=inner_cv, scoring="roc_auc",
                            max_iter=4000, random_state=42)
clf.fit(Z.values, y_25507)
P_refit = clf.predict_proba(Z.values)[:, 1]
fpr_r, tpr_r, _ = roc_curve(y_25507, P_refit)
auc_refit = g25["auc_refit_in_GSE25507"]
axB.plot(fpr_r, tpr_r, color=fs.COL["gse25507"], linewidth=1.5, linestyle=":")
axB.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axB.set_xlabel("False positive rate"); axB.set_ylabel("True positive rate")
axB.set_title("B. GSE25507 (tissue boundary)",
              loc="left", fontweight="bold", fontsize=9)
axB.set_xlim(-0.02, 1.02); axB.set_ylim(-0.02, 1.02)
axB.set_xticks([0, 0.5, 1.0]); axB.set_yticks([0, 0.5, 1.0])
h_25 = mlines.Line2D([0],[0], color=fs.COL["gse25507"], linewidth=1.5,
                       label=f"Locked panel applied (n=146, lymphocyte-enriched):\n"
                             f"AUC = {auc:.3f}  [{ci[0]:.3f}, {ci[1]:.3f}],  "
                             f"permutation p = 0.46")
h_refit = mlines.Line2D([0],[0], color=fs.COL["gse25507"], linewidth=1.5,
                          linestyle=":",
                          label=f"Within-tissue refit (same 7 genes):  "
                                f"AUC = {auc_refit:.3f} (apparent)")
axB.legend(handles=[h_25, h_refit, h_diag], loc="upper center",
            bbox_to_anchor=(0.5, -0.25), ncol=1, frameon=False,
            fontsize=6.8, handlelength=2.0)

fs.save_fig(fig, "Figure_5")
plt.close(fig)
print("Figure 5 done.")
