"""make_fig6.py v6 — Figure 6: monocyte tracking (was Fig 5 in v5)."""
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import fig_style as fs

model7 = fs.load_locked_model("../locked_model_main_7gene.json")
model6 = fs.load_locked_model("../locked_model_reduced_6gene.json")

expr_train, y_train, _ = fs.load_cohort("GPL570")
expr_6244, y_6244, _ = fs.load_cohort("GPL6244")
expr_25507, y_25507, _ = fs.load_cohort("GSE25507")
expr_42133, y_42133, _ = fs.load_cohort("GSE42133")

_, LP_train = fs.predict_strict(model7, expr_train)
_, LP_6244 = fs.predict_adapted(model7, expr_6244)
_, LP_25507 = fs.predict_adapted(model7, expr_25507)
_, LP_42133 = fs.predict_adapted(model6, expr_42133)

with open("../cell_marker_analysis.json") as f:
    cm = json.load(f)
with open("../cell_marker_GSE42133.json") as f:
    cm_42133 = json.load(f)

fig = plt.figure(figsize=(fs.WIDTH_FULL, 8.0))
outer = gridspec.GridSpec(nrows=2, ncols=1, figure=fig,
                          left=0.08, right=0.96, top=0.90, bottom=0.06,
                          height_ratios=[1.0, 1.0], hspace=0.60)

# Top: 2x2 scatter
gsA = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=2, subplot_spec=outer[0, 0],
                                        wspace=0.32, hspace=0.65)
mono_markers = fs.CELL_MARKERS["Monocytes"]
cohorts = [
    ("GPL570", "GPL570 (training)", expr_train, y_train, LP_train, fs.COL["training"]),
    ("GPL6244", "GPL6244 (external)", expr_6244, y_6244, LP_6244, fs.COL["gpl6244"]),
    ("GSE42133", "GSE42133 (independent)", expr_42133, y_42133, LP_42133, fs.COL["gse42133"]),
    ("GSE25507", "GSE25507 (tissue boundary)", expr_25507, y_25507, LP_25507, fs.COL["gse25507"]),
]
json_key_map = {
    "GPL570": "GPL570 (whole blood, training)",
    "GPL6244": "GPL6244 (whole blood, ext)",
    "GSE25507": "GSE25507 (lymphocytes)",
}
for idx, (key, title, X, y, LP, col) in enumerate(cohorts):
    ax = fig.add_subplot(gsA[idx // 2, idx % 2])
    mono = fs.cell_marker_score(X, mono_markers)
    asd_mask = y == 1
    ax.scatter(mono[~asd_mask], LP[~asd_mask], s=10, alpha=0.55,
                c=fs.COL["ctrl"], edgecolors="none")
    ax.scatter(mono[asd_mask], LP[asd_mask], s=10, alpha=0.55,
                c=fs.COL["asd"], edgecolors="none")
    z = np.polyfit(mono, LP, 1)
    xs = np.linspace(mono.min(), mono.max(), 50)
    ax.plot(xs, z[0] * xs + z[1], color=col, linewidth=1.2, alpha=0.85)
    if key == "GSE42133":
        r = cm_42133["Monocytes"]["r"]
    else:
        r = cm["correlations"][json_key_map[key]]["Monocytes"]["r"]
    ax.set_title(f"{title}    r = {r:+.2f}", fontsize=8, loc="left", pad=4)
    if idx >= 2: ax.set_xlabel("Monocyte marker score", fontsize=7.5)
    if idx % 2 == 0: ax.set_ylabel("LP$_{7}$ (logit)", fontsize=7.5)

fig.text(0.08, 0.96,
         "A. LP$_{7}$ vs monocyte-marker score — positive in all four cohorts",
         ha="left", va="top", fontsize=10, fontweight="bold")
h_ctrl = mlines.Line2D([0], [0], marker="o", color="w",
                        markerfacecolor=fs.COL["ctrl"], markersize=6,
                        linewidth=0, label="Typically developing control")
h_asd = mlines.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=fs.COL["asd"], markersize=6,
                       linewidth=0, label="ASD")
h_fit = mlines.Line2D([0], [0], color="#555", linewidth=1.2,
                       label="Linear fit (per cohort)")
fig.legend(handles=[h_ctrl, h_asd, h_fit], loc="upper center",
           bbox_to_anchor=(0.5, 0.50), ncol=3, frameon=False, fontsize=7.5,
           handlelength=1.8, handletextpad=0.5)

# Bottom: heatmap
axB = fig.add_subplot(outer[1, 0])
cohort_labels = ["GPL570\n(training,\nwhole blood)",
                 "GPL6244\n(external,\nwhole blood)",
                 "GSE42133\n(independent,\nleukocyte)",
                 "GSE25507\n(lymphocyte-\nenriched)"]
cohort_keys_in_json = ["GPL570 (whole blood, training)",
                       "GPL6244 (whole blood, ext)", "GSE42133",
                       "GSE25507 (lymphocytes)"]
mat = np.full((len(fs.CELL_TYPE_ORDER), 4), np.nan)
for i, ct in enumerate(fs.CELL_TYPE_ORDER):
    for j, k in enumerate(cohort_keys_in_json):
        if k == "GSE42133":
            info = cm_42133.get(ct)
        else:
            info = cm["correlations"].get(k, {}).get(ct)
        if info is not None: mat[i, j] = info["r"]

im = axB.imshow(mat, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
axB.set_yticks(range(len(fs.CELL_TYPE_ORDER)))
axB.set_yticklabels(fs.CELL_TYPE_ORDER, fontsize=8)
axB.set_xticks(range(4)); axB.set_xticklabels(cohort_labels, fontsize=7.5)
axB.set_title("B. Pearson correlation of LP$_{7}$ with cell-type marker scores",
              loc="left", fontweight="bold", fontsize=10, pad=8)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat[i, j]
        if np.isnan(v):
            axB.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="#aaa")
        else:
            color = "white" if abs(v) > 0.35 else "black"
            axB.text(j, i, f"{v:+.2f}", ha="center", va="center",
                      fontsize=8, color=color)
cbar = fig.colorbar(im, ax=axB, fraction=0.025, pad=0.02)
cbar.set_label("Pearson r", fontsize=8); cbar.ax.tick_params(labelsize=7)
axB.set_xticks(np.arange(-0.5, 4, 1), minor=True)
axB.set_yticks(np.arange(-0.5, len(fs.CELL_TYPE_ORDER), 1), minor=True)
axB.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
axB.tick_params(which="minor", length=0)

fs.save_fig(fig, "Figure_6")
plt.close(fig)
print("Figure 6 done.")
