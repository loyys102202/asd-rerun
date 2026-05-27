"""
make_graphical_abstract.py — Graphical abstract for iScience submission.

Visual summary of the central finding: LP7 tracks monocytes across all four
pediatric ASD blood transcriptome cohorts (r = +0.41 to +0.57), but its
ASD discrimination is platform- and cohort-contingent (AUC 0.76 internal,
0.69 within-series external, 0.42 independent external, 0.46 tissue boundary).

Output: SVG (editable) + PNG (raster) at 1100 × 800 px target size.
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

import fig_style as fs

# ============================================================================
# Load reference numbers
# ============================================================================
with open("../cell_marker_analysis.json") as f:
    cm = json.load(f)
with open("../cell_marker_GSE42133.json") as f:
    cm_42133 = json.load(f)

# Monocyte r across 4 cohorts
mono_r = {
    "GPL570":   cm["correlations"]["GPL570 (whole blood, training)"]["Monocytes"]["r"],
    "GPL6244":  cm["correlations"]["GPL6244 (whole blood, ext)"]["Monocytes"]["r"],
    "GSE42133": cm_42133["Monocytes"]["r"],
    "GSE25507": cm["correlations"]["GSE25507 (lymphocytes)"]["Monocytes"]["r"],
}

# AUCs across 4 cohorts (manuscript reference numbers)
aucs = {
    "GPL570":   0.762,    # internal nested CV
    "GPL6244":  0.689,
    "GSE42133": 0.420,
    "GSE25507": 0.464,
}

# ============================================================================
# Figure — designed for ~1100 × 800 px at 110 dpi → 10 × 7.3 in
# ============================================================================
fig = plt.figure(figsize=(10.0, 7.0))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# ---------- Title ----------
ax.text(50, 96,
        "A monocyte-tracking blood expression score reveals platform-specific",
        ha="center", va="top", fontsize=14, fontweight="bold", color="#222")
ax.text(50, 91,
        "dissociation between cellular composition and autism diagnosis",
        ha="center", va="top", fontsize=14, fontweight="bold", color="#222")

# ---------- 4 cohort schematic row ----------
def cohort_pill(x, y, w, h, name, n, role, role_col, body_col="#f8f9fa"):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=1.5",
                          linewidth=1.0, edgecolor=role_col, facecolor=body_col)
    ax.add_patch(box)
    # Top role strip
    badge = FancyBboxPatch((x + 0.4, y + h - 4), w - 0.8, 4,
                            boxstyle="round,pad=0.0,rounding_size=0.8",
                            linewidth=0, facecolor=role_col, alpha=0.92)
    ax.add_patch(badge)
    ax.text(x + w / 2, y + h - 2, role, ha="center", va="center",
            color="white", fontsize=8, fontweight="bold")
    ax.text(x + w / 2, y + h - 7.5, name, ha="center", va="center",
            fontsize=8.5, fontweight="bold", color="#333")
    ax.text(x + w / 2, y + h - 10.5, n, ha="center", va="center",
            fontsize=8, color="#555")

cohorts_top = [
    (3.5, 75, 21.5, 13, "GSE18123-GPL570", "n = 99\nwhole blood",
     "Training", fs.COL["training"]),
    (27, 75, 21.5, 13, "GSE18123-GPL6244", "n = 186\nwhole blood",
     "External", fs.COL["gpl6244"]),
    (51, 75, 21.5, 13, "GSE42133", "n = 147\nleukocyte",
     "Independent", fs.COL["gse42133"]),
    (75, 75, 21.5, 13, "GSE25507", "n = 146\nlymphocyte",
     "Tissue boundary", fs.COL["gse25507"]),
]
for c in cohorts_top:
    cohort_pill(*c)

# ---------- Locked panel pill ----------
panel_y = 64
panel_box = FancyBboxPatch((15, panel_y - 0.5), 70, 6.0,
                            boxstyle="round,pad=0.0,rounding_size=1.5",
                            linewidth=1.0, edgecolor="#444",
                            facecolor="#eef4fa")
ax.add_patch(panel_box)
ax.text(50, panel_y + 2.5,
        "Locked 7-gene panel  •  stability-selected LASSO  •  "
        "KYNU  CMYA5  MAPK8IP1  CES1  KIAA0087  POU2AF1  ERVK3-2",
        ha="center", va="center", fontsize=8.5, color="#333", fontweight="bold")

# Arrow from cohort row → panel
for cx, cy, cw, ch, *_ in cohorts_top:
    arr = FancyArrowPatch((cx + cw / 2, cy), (cx + cw / 2 + (50 - cx - cw / 2) * 0.0, panel_y + 5.5),
                          arrowstyle="->,head_width=2.5,head_length=2.5",
                          linewidth=0.8, color="#888")
    ax.add_patch(arr)

# ---------- Two-track outcome row ----------
# Headers
ax.text(25, 55,
        "What it MEASURES",
        ha="center", va="center", fontsize=10, fontweight="bold",
        color="#0b6e4f")
ax.text(75, 55,
        "What it DISCRIMINATES",
        ha="center", va="center", fontsize=10, fontweight="bold",
        color="#a52a2a")

# Left track: monocyte tracking bar chart
left_x = 7
right_x = 43
y_top = 50
y_bot = 22
cohort_names = ["GPL570", "GPL6244", "GSE42133", "GSE25507"]
cohort_labels_short = ["GPL570\n(training)", "GPL6244\n(external)",
                       "GSE42133\n(indep.)", "GSE25507\n(tissue)"]
colors_per_cohort = [fs.COL["training"], fs.COL["gpl6244"],
                      fs.COL["gse42133"], fs.COL["gse25507"]]

# Plot the monocyte r values as a horizontal mini-chart
# Use a sub-axes overlay via inset_axes for cleaner rendering
ax_mono = fig.add_axes([0.07, 0.30, 0.36, 0.20])  # [left, bottom, width, height]
ax_mono.set_facecolor("none")
mono_vals = [mono_r[k] for k in cohort_names]
ypos = np.arange(len(cohort_names))[::-1]
bars = ax_mono.barh(ypos, mono_vals, color=colors_per_cohort,
                      edgecolor="black", linewidth=0.6, height=0.65)
ax_mono.axvline(0, color="black", linewidth=0.6)
ax_mono.set_yticks(ypos)
ax_mono.set_yticklabels(cohort_labels_short, fontsize=8)
ax_mono.set_xlim(-0.05, 0.75)
ax_mono.set_xlabel("Pearson r  (LP$_7$ vs monocyte score)", fontsize=8.5)
for spine in ["top", "right"]:
    ax_mono.spines[spine].set_visible(False)
for bar, v in zip(bars, mono_vals):
    ax_mono.text(v + 0.02, bar.get_y() + bar.get_height() / 2,
                 f"+{v:.2f}", va="center", ha="left", fontsize=8.5,
                 fontweight="bold")
ax_mono.set_title("Monocyte signal preserved across all 4 cohorts",
                   loc="center", fontsize=8.5, color="#0b6e4f", pad=4)

# Right track: AUC bars
ax_auc = fig.add_axes([0.55, 0.30, 0.36, 0.20])
ax_auc.set_facecolor("none")
auc_vals = [aucs[k] for k in cohort_names]
ypos = np.arange(len(cohort_names))[::-1]
bars = ax_auc.barh(ypos, auc_vals, color=colors_per_cohort,
                    edgecolor="black", linewidth=0.6, height=0.65)
ax_auc.axvline(0.5, color="#888", linewidth=0.8, linestyle="--",
                label="Chance")
ax_auc.set_yticks(ypos)
ax_auc.set_yticklabels(cohort_labels_short, fontsize=8)
ax_auc.set_xlim(0.0, 1.0)
ax_auc.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
ax_auc.set_xlabel("AUC for ASD discrimination", fontsize=8.5)
for spine in ["top", "right"]:
    ax_auc.spines[spine].set_visible(False)
for bar, v in zip(bars, auc_vals):
    # Annotate inside bar end (always to the right of bar)
    ax_auc.text(v + 0.01, bar.get_y() + bar.get_height() / 2,
                 f"{v:.2f}", va="center", ha="left", fontsize=8.5,
                 fontweight="bold")
ax_auc.set_title("ASD discrimination is platform-contingent",
                  loc="center", fontsize=8.5, color="#a52a2a", pad=4)

# ---------- Conclusion strip ----------
concl_y = 7
concl_box = FancyBboxPatch((4, concl_y - 1), 92, 8.5,
                            boxstyle="round,pad=0.0,rounding_size=1.5",
                            linewidth=1.2, edgecolor="#444",
                            facecolor="#fff5e6")
ax.add_patch(concl_box)
ax.text(50, concl_y + 5.5,
        "The blood-expression signature primarily measures monocyte-compartment biology,",
        ha="center", va="center", fontsize=10, fontweight="bold", color="#333")
ax.text(50, concl_y + 2.5,
        "while its association with ASD is contingent on cohort- and platform-specific "
        "monocyte–diagnosis coupling.",
        ha="center", va="center", fontsize=10, fontweight="bold", color="#333")

# Save SVG (editable) + PNG (raster) + PDF
from pathlib import Path
out_dir = Path("output")
out_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(out_dir / "Graphical_Abstract.svg", format="svg")
fig.savefig(out_dir / "Graphical_Abstract.png", dpi=300)
fig.savefig(out_dir / "Graphical_Abstract.pdf")
plt.close(fig)
print("Graphical abstract done — wrote SVG, PNG (300 dpi), PDF.")
