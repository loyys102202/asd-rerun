"""
make_fig4.py v6 — Figure 4: external validation on GSE18123-GPL6244.

(Was Figure 3 in the previous 5-figure layout.)
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from sklearn.metrics import roc_curve, roc_auc_score, brier_score_loss

import fig_style as fs

model = fs.load_locked_model("../locked_model_main_7gene.json")
expr, y, pheno = fs.load_cohort("GPL6244")
P_adapt, LP_adapt = fs.predict_adapted(model, expr)
P_strict, LP_strict = fs.predict_strict(model, expr)

with open("../calibration.json") as f:
    cal = json.load(f)
alpha_recal = cal["external_adapted"]["cal_intercept"]
beta_recal = cal["external_adapted"]["cal_slope"]

with open("../robustness_checks.json") as f:
    rc = json.load(f)

fig = plt.figure(figsize=(fs.WIDTH_FULL, 6.5))
gs = gridspec.GridSpec(
    nrows=2, ncols=2,
    left=0.07, right=0.97, top=0.94, bottom=0.07,
    wspace=0.40, hspace=0.80,
)

# Panel A
axA = fig.add_subplot(gs[0, 0])
fpr_a, tpr_a, _ = roc_curve(y, P_adapt)
fpr_s, tpr_s, _ = roc_curve(y, P_strict)
auc_a = roc_auc_score(y, P_adapt); auc_s = roc_auc_score(y, P_strict)
ci_a = fs.delong_ci(y, P_adapt);   ci_s = fs.delong_ci(y, P_strict)

axA.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axA.plot(fpr_a, tpr_a, color=fs.COL["adapted"], linewidth=1.5)
axA.plot(fpr_s, tpr_s, color=fs.COL["strict"], linewidth=1.5, linestyle=":")
axA.set_xlabel("False positive rate"); axA.set_ylabel("True positive rate")
axA.set_title("A. GSE18123-GPL6244 ROC", loc="left", fontweight="bold", fontsize=9)
axA.set_xlim(-0.02, 1.02); axA.set_ylim(-0.02, 1.02)
axA.set_xticks([0, 0.5, 1.0]); axA.set_yticks([0, 0.5, 1.0])
h_a = mlines.Line2D([0],[0], color=fs.COL["adapted"], linewidth=1.5,
                     label=f"Platform-adapted: AUC = {auc_a:.3f}  [{ci_a[0]:.3f}, {ci_a[1]:.3f}]")
h_s = mlines.Line2D([0],[0], color=fs.COL["strict"], linewidth=1.5, linestyle=":",
                     label=f"Strict train-scale: AUC = {auc_s:.3f}  [{ci_s[0]:.3f}, {ci_s[1]:.3f}]")
axA.legend(handles=[h_a, h_s], loc="upper center", bbox_to_anchor=(0.5, -0.22),
            ncol=1, frameon=False, fontsize=6.8, handlelength=2.0)

# Panel B
axB = fig.add_subplot(gs[0, 1])
pred_pre, obs_pre, _ = fs.calibration_bin(y, P_adapt, n_bins=10, strategy="quantile")
LP_recal = alpha_recal + beta_recal * LP_adapt
P_recal = 1.0 / (1.0 + np.exp(-np.clip(LP_recal, -50, 50)))
pred_post, obs_post, _ = fs.calibration_bin(y, P_recal, n_bins=10, strategy="quantile")
brier_pre = brier_score_loss(y, P_adapt); brier_post = brier_score_loss(y, P_recal)

axB.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axB.plot(pred_pre, obs_pre, "o-", color=fs.COL["before"], linewidth=1.4,
         markersize=4.5, zorder=2)
axB.plot(pred_post, obs_post, "s-", color=fs.COL["after"], linewidth=1.4,
         markersize=4.5, zorder=3)
axB.set_xlabel("Predicted probability of ASD"); axB.set_ylabel("Observed proportion ASD")
axB.set_title("B. Calibration before / after recalibration",
              loc="left", fontweight="bold", fontsize=9)
axB.set_xlim(-0.02, 1.02); axB.set_ylim(-0.02, 1.02)
axB.set_xticks([0, 0.5, 1.0]); axB.set_yticks([0, 0.5, 1.0])
h_pre = mlines.Line2D([0],[0], color=fs.COL["before"], marker="o", linewidth=1.4,
                       markersize=4, label=f"Before  (Brier = {brier_pre:.3f})")
h_post = mlines.Line2D([0],[0], color=fs.COL["after"], marker="s", linewidth=1.4,
                        markersize=4, label=f"After   (Brier = {brier_post:.3f})")
h_ideal = mlines.Line2D([0],[0], color=fs.COL["diag"], linestyle="--",
                          linewidth=0.7, label="Ideal")
axB.legend(handles=[h_pre, h_post, h_ideal], loc="upper center",
            bbox_to_anchor=(0.5, -0.22), ncol=1, frameon=False,
            fontsize=6.8, handlelength=2.0)
axB.annotate(f"Recalibration: α = {alpha_recal:+.3f},  β = {beta_recal:+.3f}",
              xy=(0.5, -0.45), xycoords="axes fraction",
              ha="center", va="top", fontsize=6.8, color="#555", style="italic")

# Panel C
axC = fig.add_subplot(gs[1, 0])
ph_idx = pheno.set_index("sample_id").loc[expr.index]
mask_m = ph_idx["sex"].astype(str).str.lower().str.startswith("m").values
mask_f = ph_idx["sex"].astype(str).str.lower().str.startswith("f").values
legend_items = []
for mask, label_name, col in [(mask_m, "Males", fs.COL["male"]),
                                (mask_f, "Females", fs.COL["female"])]:
    if mask.sum() < 5 or len(np.unique(y[mask])) < 2: continue
    P_sub = P_adapt[mask]; y_sub = y[mask]
    fpr, tpr, _ = roc_curve(y_sub, P_sub)
    auc_sub = roc_auc_score(y_sub, P_sub); ci_sub = fs.delong_ci(y_sub, P_sub)
    n_asd = int((y_sub == 1).sum()); n_ctl = int((y_sub == 0).sum())
    axC.plot(fpr, tpr, color=col, linewidth=1.5)
    legend_items.append(mlines.Line2D(
        [0], [0], color=col, linewidth=1.5,
        label=f"{label_name} (n={mask.sum()}; {n_asd} ASD / {n_ctl} ctrl):  "
              f"AUC = {auc_sub:.3f}  [{ci_sub[0]:.3f}, {ci_sub[1]:.3f}]"))
axC.plot([0, 1], [0, 1], "--", color=fs.COL["diag"], linewidth=0.7)
axC.set_xlabel("False positive rate"); axC.set_ylabel("True positive rate")
axC.set_title("C. Sex-stratified ROC (platform-adapted)",
              loc="left", fontweight="bold", fontsize=9)
axC.set_xlim(-0.02, 1.02); axC.set_ylim(-0.02, 1.02)
axC.set_xticks([0, 0.5, 1.0]); axC.set_yticks([0, 0.5, 1.0])
axC.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.22),
            ncol=1, frameon=False, fontsize=6.8, handlelength=2.0)

# Panel D
axD = fig.add_subplot(gs[1, 1])
null_mean = rc["perm_external_null_mean"]
null_lo, null_hi = rc["perm_external_null_95"]
observed_auc = rc["observed_external_auc"]
null_sd = (null_hi - null_lo) / (2 * 1.96)
xs = np.linspace(0.35, 0.85, 400)
ys = (1.0 / (null_sd * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((xs - null_mean) / null_sd) ** 2)
ys = ys / ys.max() * 1.0
axD.fill_between(xs, ys, 0, color=fs.COL["null"], alpha=0.45)
axD.fill_between(xs[(xs >= null_lo) & (xs <= null_hi)],
                 ys[(xs >= null_lo) & (xs <= null_hi)], 0,
                 color=fs.COL["null"], alpha=0.30)
axD.axvline(null_lo, color="#5a5a5a", linewidth=0.6, linestyle=":")
axD.axvline(null_hi, color="#5a5a5a", linewidth=0.6, linestyle=":")
axD.axvline(observed_auc, color=fs.COL["obs"], linewidth=1.6)
axD.set_xlim(0.35, 0.85); axD.set_ylim(0, 1.15)
axD.set_xlabel("AUC under null"); axD.set_ylabel("Density (normalised)")
axD.set_title("D. Permutation null", loc="left", fontweight="bold", fontsize=9)
h_null = mpatches.Patch(facecolor=fs.COL["null"], alpha=0.45,
                         label=f"Permutation null (10,000 perms; mean = {null_mean:.3f})")
h_band = mpatches.Patch(facecolor=fs.COL["null"], alpha=0.30,
                         label=f"95% null band  [{null_lo:.3f}, {null_hi:.3f}]")
h_obs = mlines.Line2D([0],[0], color=fs.COL["obs"], linewidth=1.6,
                       label=f"Observed AUC = {observed_auc:.3f}  (p < 0.0001)")
axD.legend(handles=[h_null, h_band, h_obs], loc="upper center",
           bbox_to_anchor=(0.5, -0.22), ncol=1, frameon=False,
           fontsize=6.8, handlelength=2.0)

fs.save_fig(fig, "Figure_4")
plt.close(fig)
print("Figure 4 done.")
