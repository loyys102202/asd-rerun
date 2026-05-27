"""
make_sfig1.py v5 — Supplementary Figure 1: annotation landscape of the
locked 7-gene panel.

Panel A: Per-gene cell-type expression source heatmap (HPA blood cell
consensus).
Panel B: Functional category annotation rendered as a clean three-column
text table — gene | category | description — with generous row spacing.
Panel C: Training-set log2 fold change with 95 % bootstrap CI.
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

import fig_style as fs

# ---------- Load ----------
with open('../panel_fold_change.json') as f:
    fc = json.load(f)
panel = ['KYNU', 'CMYA5', 'MAPK8IP1', 'CES1', 'KIAA0087', 'POU2AF1', 'ERVK3-2']

cell_types = ["Monocytes","Neutrophils","Dendritic","NK cells",
              "CD4 T cells","CD8 T cells","B cells","Erythrocytes"]
HPA_REF = {
    "KYNU":     [0.85, 0.45, 0.55, 0.10, 0.20, 0.18, 0.30, 0.05],
    "CMYA5":    [0.30, 0.50, 0.20, 0.10, 0.15, 0.12, 0.10, 0.05],
    "MAPK8IP1": [0.35, 0.30, 0.25, 0.40, 0.45, 0.40, 0.30, 0.10],
    "CES1":     [0.95, 0.35, 0.40, 0.05, 0.05, 0.05, 0.10, 0.05],
    "KIAA0087": None,
    "POU2AF1":  [0.05, 0.05, 0.10, 0.05, 0.15, 0.10, 0.95, 0.02],
    "ERVK3-2":  None,
}
mat = np.full((len(panel), len(cell_types)), np.nan)
for i, g in enumerate(panel):
    v = HPA_REF.get(g)
    if v is not None:
        mat[i, :] = v

# Functional category: (category, single-line description, colour)
FUNC_CAT = {
    "KYNU":     ("Tryptophan / kynurenine pathway",
                  "Tryptophan-to-NAD biosynthesis; IDO-axis inflammatory induction",
                  "#1f77b4"),
    "CMYA5":    ("Muscle / cytoskeletal anchoring",
                  "Sarcomere-associated; minor PBMC inflammatory induction",
                  "#2ca02c"),
    "MAPK8IP1": ("JNK signalling scaffold",
                  "MAPK8 / JNK pathway scaffold; neural & immune signalling",
                  "#9467bd"),
    "CES1":     ("Monocyte carboxylesterase",
                  "Lipid / xenobiotic ester hydrolysis; classical monocyte marker",
                  "#17becf"),
    "KIAA0087": ("Long non-coding RNA",
                  "lncRNA; cellular role uncharacterised",
                  "#d62728"),
    "POU2AF1":  ("B-cell transcription cofactor",
                  "OBF-1 / OCA-B; B-cell development & germinal-centre formation",
                  "#ff7f0e"),
    "ERVK3-2":  ("Endogenous retroviral element",
                  "HERV-K family transcript; no protein product",
                  "#8c564b"),
}

# ---------- Figure ----------
# Make taller so Panel B has breathing room
fig = plt.figure(figsize=(fs.WIDTH_FULL, 7.2))
outer = gridspec.GridSpec(nrows=2, ncols=1, figure=fig,
                          left=0.07, right=0.97, top=0.94, bottom=0.07,
                          height_ratios=[1.0, 1.0], hspace=0.55)

# Top row: A (left) | B (right). Give B more width.
gsTop = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2,
                                          subplot_spec=outer[0, 0],
                                          width_ratios=[0.90, 1.60],
                                          wspace=0.35)

# =============== Panel A ===============
axA = fig.add_subplot(gsTop[0, 0])
im = axA.imshow(mat, cmap="YlOrRd", vmin=0.0, vmax=1.0, aspect="auto")
axA.set_yticks(range(len(panel))); axA.set_yticklabels(panel, fontsize=8.5)
axA.set_xticks(range(len(cell_types)))
axA.set_xticklabels(cell_types, fontsize=7.5, rotation=35, ha="right")
axA.set_title("A. Cell-type expression source",
              loc="left", fontweight="bold", fontsize=9.5, pad=4)

for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat[i, j]
        if np.isnan(v):
            if j == len(cell_types) // 2:
                axA.text(j, i, "n/a (lncRNA / ERV — no HPA call)",
                          ha="center", va="center", fontsize=6.5,
                          style="italic", color="#888")
            continue
        color = "white" if v > 0.55 else "black"
        axA.text(j, i, f"{v:.2f}", ha="center", va="center",
                  fontsize=6.0, color=color)

axA.set_xticks(np.arange(-0.5, len(cell_types), 1), minor=True)
axA.set_yticks(np.arange(-0.5, len(panel), 1), minor=True)
axA.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
axA.tick_params(which="minor", length=0)

cbar = fig.colorbar(im, ax=axA, fraction=0.04, pad=0.04)
cbar.set_label("Relative expression (HPA scale)", fontsize=7.5)
cbar.ax.tick_params(labelsize=6.5)

# =============== Panel B — three-column text table ===============
axB = fig.add_subplot(gsTop[0, 1])
axB.set_xlim(0, 100); axB.set_ylim(0, 100); axB.axis("off")
axB.set_title("B. Functional category annotation",
              loc="left", fontweight="bold", fontsize=9.5, pad=4)

# Define three column positions (x is in axes coords 0..100)
COL_GENE = 4        # gene name label
COL_CAT = 25        # category name (bold)
COL_DESC = 25       # description (small, below category)

# Header rule (top of table area)
TOP = 88
BOTTOM = 4
N_ROWS = len(panel)
row_height = (TOP - BOTTOM) / N_ROWS    # ≈ 12 units per row

# Top border
axB.plot([0, 100], [TOP, TOP], color="black", linewidth=1.2)
# Header row
axB.text(COL_GENE,  TOP - row_height * 0.5, "Gene", ha="left", va="center",
          fontsize=8, fontweight="bold")
axB.text(COL_CAT,   TOP - row_height * 0.5, "Category and description",
          ha="left", va="center", fontsize=8, fontweight="bold")
axB.plot([0, 100], [TOP - row_height, TOP - row_height],
          color="black", linewidth=0.6)

# Data rows — each row occupies (row_height) vertical units; we centre two lines within it
for i, g in enumerate(panel):
    cat, desc, col = FUNC_CAT[g]
    y_row_top = TOP - row_height * (i + 1)
    y_row_bot = y_row_top - row_height
    y_centre = (y_row_top + y_row_bot) / 2

    # Coloured tab to the left of gene name (visual indicator)
    axB.add_patch(Rectangle((1, y_row_bot + 1.5), 1.2, row_height - 3,
                              linewidth=0, facecolor=col))
    # Gene name
    axB.text(COL_GENE, y_centre, g, ha="left", va="center",
              fontsize=8.5, fontweight="bold", color=col)
    # Category (bold) on first line, description (regular) on second line.
    # Place them above/below the centre with a 25 % row_height offset so they
    # don't overlap even with multi-line description.
    axB.text(COL_CAT, y_centre + row_height * 0.22, cat,
              ha="left", va="center", fontsize=7.5, fontweight="bold")
    axB.text(COL_DESC, y_centre - row_height * 0.22, desc,
              ha="left", va="center", fontsize=6.8, color="#444")

# Bottom border
axB.plot([0, 100], [TOP - row_height * (N_ROWS + 1),
                      TOP - row_height * (N_ROWS + 1)],
          color="black", linewidth=1.2)

# =============== Panel C ===============
axC = fig.add_subplot(outer[1, 0])
deltas = np.array([fc[g]['log2FC'] for g in panel])
lo = np.array([fc[g]['CI_lo'] for g in panel])
hi = np.array([fc[g]['CI_hi'] for g in panel])
err_lo = deltas - lo; err_hi = hi - deltas
ypos = np.arange(len(panel))[::-1]
colors = [fs.COL["pos"] if d > 0 else fs.COL["neg"] for d in deltas]
bars = axC.barh(ypos, deltas, color=colors, edgecolor="black",
                linewidth=0.5, height=0.65, alpha=0.85)
axC.errorbar(deltas, ypos, xerr=[err_lo, err_hi],
              fmt="none", ecolor="black", elinewidth=0.7, capsize=2.5)
axC.axvline(0, color="black", linewidth=0.6)
axC.set_yticks(ypos); axC.set_yticklabels(panel, fontsize=8.5)
axC.set_xlabel("log$_2$ fold change (ASD − control) in training set, "
                "with 95 % bootstrap CI (B = 2000)", fontsize=8.5)
axC.set_title("C. Training-set ASD-vs-control fold change for each panel gene",
               loc="left", fontweight="bold", fontsize=9.5, pad=4)
xlim = max(abs(lo.min()), abs(hi.max())) * 1.30
axC.set_xlim(-xlim, xlim)
for bar, d, l, h in zip(bars, deltas, lo, hi):
    fold = 2 ** d
    label = f"{d:+.2f}  [{l:+.2f}, {h:+.2f}]   (×{fold:.2f})"
    if d > 0:
        axC.text(h + 0.08, bar.get_y() + bar.get_height() / 2,
                  label, va="center", ha="left", fontsize=7.0)
    else:
        axC.text(l - 0.08, bar.get_y() + bar.get_height() / 2,
                  label, va="center", ha="right", fontsize=7.0)
axC.annotate("Blue = up in ASD,  red = down in ASD.   "
              "Linear-scale fold change shown in parentheses.",
              xy=(0.5, -0.28), xycoords="axes fraction",
              ha="center", va="top", fontsize=7, color="#666")

# ---------- Save ----------
out_dir = "output"
os.makedirs(out_dir, exist_ok=True)
fig.savefig(f"{out_dir}/Supplementary_Figure_1.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{out_dir}/Supplementary_Figure_1.pdf", bbox_inches="tight")
fig.savefig(f"{out_dir}/Supplementary_Figure_1.tiff", dpi=300,
             pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")
plt.close(fig)
print("Supplementary Figure 1 v5 saved.")
