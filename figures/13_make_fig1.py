"""
make_fig1.py v7 — Figure 1: study design (standalone).

Arrows route below the cohort row with small curvature so they stay well
above the bottom caption.  Boxes widened to fit role labels.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

import fig_style as fs

fig = plt.figure(figsize=(fs.WIDTH_FULL, 4.2))
ax = fig.add_axes([0.04, 0.04, 0.94, 0.88])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

ax.text(0, 99, "Figure 1.  Study design — four pediatric ASD blood transcriptome cohorts",
        ha="left", va="top", fontsize=11, fontweight="bold")

def cohort_box(x, y, w, h, lines, color, label_role):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.0,rounding_size=1.5",
                         linewidth=1.4, edgecolor=color, facecolor="white",
                         zorder=3)
    ax.add_patch(box)
    badge = FancyBboxPatch((x + 0.5, y + h - 7), w - 1, 7,
                           boxstyle="round,pad=0.0,rounding_size=0.8",
                           linewidth=0, facecolor=color, alpha=0.92, zorder=4)
    ax.add_patch(badge)
    ax.text(x + w / 2, y + h - 3.5, label_role, ha="center", va="center",
            color="white", fontsize=7.5, fontweight="bold", zorder=5)
    for i, line in enumerate(lines):
        ax.text(x + w / 2, y + h - 14 - i * 5.5, line, ha="center", va="center",
                fontsize=7.8, zorder=5)

# Wider boxes (23 each), 4 boxes fit width 0..100 with 1.5 padding between
xs = [2, 25.5, 49, 72.5]
yb, w, h = 32, 23, 56
cohorts = [
    (xs[0], yb, w, h,
     ["GSE18123-GPL570","Affymetrix U133+2","Whole blood",
      "n = 99","66 ASD / 33 ctrl","All male"],
     fs.COL["training"], "Training"),
    (xs[1], yb, w, h,
     ["GSE18123-GPL6244","Affy HuGene 1.0 ST","Whole blood",
      "n = 186","104 ASD / 82 ctrl","128 M / 58 F"],
     fs.COL["gpl6244"], "External (same series)"),
    (xs[2], yb, w, h,
     ["GSE42133","Illumina HT-12 v4","Leukocyte",
      "n = 147","91 ASD / 56 ctrl","All male"],
     fs.COL["gse42133"], "Independent external"),
    (xs[3], yb, w, h,
     ["GSE25507","Affymetrix U133+2","Lymphocyte-enriched",
      "n = 146","82 ASD / 64 ctrl","Sex n/a"],
     fs.COL["gse25507"], "Tissue boundary"),
]
for c in cohorts:
    cohort_box(*c)

# Arrows below boxes — gentle arc so they stay well above caption
training_bot = (xs[0] + w/2, yb)
for x_target in xs[1:]:
    target_bot = (x_target + w/2, yb)
    arr = FancyArrowPatch(training_bot, target_bot,
                          connectionstyle="arc3,rad=0.20",
                          arrowstyle="->,head_width=2.8,head_length=3.4",
                          linewidth=1.2, color="#777", zorder=2)
    ax.add_patch(arr)

# Caption strip — explicit two-line format
ax.text(50, 8,
        "Locked panel learnt on GSE18123-GPL570, then applied to the three external / boundary cohorts.   Total n = 578.",
        ha="center", va="center", fontsize=7.8, style="italic", color="#555")
ax.text(50, 3,
        "Unified preprocessing: log$_2$ → IQR-max probe-to-gene → GPL570 ∩ GPL6244 common-gene pool (17,923 genes).",
        ha="center", va="center", fontsize=7.8, style="italic", color="#555")

fs.save_fig(fig, "Figure_1")
plt.close(fig)
print("Figure 1 done.")
