"""
make_fig2.py v6 — Figure 2: locked 7-gene panel.

Standalone figure for the locked panel: coefficients (left) + stability
selection frequencies (right).
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig_style as fs

model = fs.load_locked_model("../locked_model_main_7gene.json")
freq = pd.read_csv("../nestedcv_gene_freq.csv", index_col=0)

panel = model["panel"]
coefs = pd.Series(model["coefficients"]).reindex(panel)
order = coefs.abs().sort_values(ascending=False).index.tolist()
coefs = coefs.reindex(order)
freq_map = freq["freq"].to_dict()
sel_freq = pd.Series({g: freq_map.get(g, np.nan) for g in order})

fig = plt.figure(figsize=(fs.WIDTH_FULL, 3.6))
gs = gridspec.GridSpec(nrows=1, ncols=2,
                       wspace=0.40,
                       left=0.07, right=0.97, top=0.90, bottom=0.20)

# ===== Panel A: coefficients =====
axA = fig.add_subplot(gs[0, 0])
ypos = np.arange(len(coefs))[::-1]
colors = [fs.COL["pos"] if v > 0 else fs.COL["neg"] for v in coefs.values]
bars = axA.barh(ypos, coefs.values, color=colors, edgecolor="black",
                linewidth=0.5, height=0.7)
axA.axvline(0, color="black", linewidth=0.6)
axA.set_yticks(ypos); axA.set_yticklabels(coefs.index, fontsize=8.5)
axA.set_xlabel("Logit-scale coefficient", fontsize=8.5)
axA.set_title("A. Locked 7-gene panel coefficients",
              loc="left", fontweight="bold", fontsize=9.5)
axA.set_xlim(-1.3, 2.0)
for bar, val in zip(bars, coefs.values):
    if val > 0:
        axA.text(val + 0.08, bar.get_y() + bar.get_height()/2,
                 f"+{val:.3f}", va="center", ha="left", fontsize=8)
    else:
        axA.text(val - 0.08, bar.get_y() + bar.get_height()/2,
                 f"{val:.3f}", va="center", ha="right", fontsize=8)
axA.text(0.97, 0.03, f"Intercept β₀ = +{model['intercept']:.4f}",
         transform=axA.transAxes, ha="right", va="bottom",
         fontsize=7.8, style="italic")
axA.annotate("Bars: blue = up in ASD,  red = down in ASD",
             xy=(0.5, -0.30), xycoords="axes fraction",
             ha="center", va="top", fontsize=7.5, color="#666")

# ===== Panel B: selection frequencies =====
axB = fig.add_subplot(gs[0, 1])
ypos = np.arange(len(sel_freq))[::-1]
bars = axB.barh(ypos, sel_freq.values * 100, color="#7f7f7f",
                edgecolor="black", linewidth=0.5, height=0.7)
axB.axvline(60, color=fs.COL["obs"], linewidth=0.9, linestyle="--")
axB.set_yticks(ypos); axB.set_yticklabels(sel_freq.index, fontsize=8.5)
axB.set_xlabel("Selection frequency across 500 outer folds (%)", fontsize=8.5)
axB.set_xlim(0, 110)
axB.set_title("B. Stability selection across 500 outer folds",
              loc="left", fontweight="bold", fontsize=9.5)
for bar, val in zip(bars, sel_freq.values):
    axB.text(val * 100 + 1.5, bar.get_y() + bar.get_height()/2,
             f"{val * 100:.1f}%", va="center", ha="left", fontsize=8)
axB.annotate("Dashed red line: 60% stability threshold",
             xy=(0.5, -0.30), xycoords="axes fraction",
             ha="center", va="top", fontsize=7.5, color="#666")

fs.save_fig(fig, "Figure_2")
plt.close(fig)
print("Figure 2 done.")
