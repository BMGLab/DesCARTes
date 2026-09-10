#!/usr/bin/env python3
"""
Composite Figure 2 (panels A-E).

  A  De novo design and off-target screening workflow
  B  STRIVE specificity landscape (regenerated as vector from the deposited scoring
     table, replacing the raster panel from scfv_specificity_decision_notebook.ipynb)
  C  Predicted scFv-73-0/CLDN4 complex
  D  VH-CLDN4 interface contacts
  E  VL-CLDN4 interface contacts

Panels A and C-E are the authors' renders, re-extracted at native resolution from the
figure sources rather than re-screenshotted. Panels C-E are raster at 144-200 ppi in the
supplied files; at print size that is below 300 dpi and they should be re-rendered from
PyMOL before final submission.

Outputs: figures/Figure2_composite.{pdf,png}
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from figcheck import assert_no_text_overlap

O, F = "output/", "../figures/"
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
TIERS = [("Not specific (<= 1.0)", "#b7d3f6"), ("Acceptable (1.0-1.2)", "#6da7ec"),
         ("Good (1.2-1.5)", "#2a78d6"), ("Excellent (> 1.5)", "#104281")]
plt.rcParams.update({"font.size":7.5,"axes.edgecolor":"#c9c8c4","axes.labelcolor":INK2,
                     "xtick.color":INK2,"ytick.color":INK2,"figure.facecolor":SURF,
                     "axes.facecolor":SURF,"savefig.facecolor":SURF})
def panel(ax, letter, dx=-0.06):
    ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="bottom", ha="left")
def imgpanel(ax, path, letter, title=None, dx=-0.04):
    ax.imshow(mpimg.imread(F + "source_figure2/" + path))
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    if title: ax.set_title("   " + title, fontsize=7.5, color=INK, loc="left")
    panel(ax, letter, dx)

fig = plt.figure(figsize=(7.4, 9.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.42, 1.28, 0.72], hspace=0.30, wspace=0.22)

# ------------------------------------------------------------------- A ------
axA = fig.add_subplot(gs[0, :]); imgpanel(axA, "A_workflow.png", "A", dx=-0.022)

# ------------------------------------------------------------------- B ------
axB = fig.add_subplot(gs[1, 0])
d = pd.read_csv(O + "TableS_specificity_full.csv")
d["label"] = d.scFv.str.replace("scfv_", "scFv-", regex=False).str.replace("_", "-", regex=False)
d["tier"] = d.R_spec.map(lambda r: 3 if r > 1.5 else 2 if r > 1.2 else 1 if r > 1.0 else 0)
axB.axhline(1.0, color="#9a9994", lw=0.8, ls="--")
axB.scatter(d.CLDN4, d.R_spec, s=26, c=d.tier.map(lambda t: TIERS[t][1]),
            edgecolor=SURF, linewidth=0.5, zorder=3)
for _, r in d[d.R_spec > 1.2].iterrows():
    axB.annotate(r.label, (r.CLDN4, r.R_spec), textcoords="offset points",
                 xytext=(-5, 5), ha="right", fontsize=6, color=INK)
axB.set_xlabel("On-target CLDN4 STRIVE score")
axB.set_ylabel("R$_{spec}$  (specificity ratio)")
axB.legend(handles=[Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor=c,
                           markeredgecolor=SURF, label=l) for l, c in TIERS],
           fontsize=5.8, frameon=False, loc="upper left", handletextpad=0.4,
           borderpad=0.15, labelspacing=0.3)
for s in ("top", "right"): axB.spines[s].set_visible(False)
panel(axB, "B", dx=-0.16)

# ------------------------------------------------------------------- C ------
axC = fig.add_subplot(gs[1, 1]); imgpanel(axC, "C_complex.png", "C", dx=-0.02)

# ---------------------------------------------------------------- D and E ---
axD = fig.add_subplot(gs[2, 0]); imgpanel(axD, "D_heavy.png", "D", "VH-CLDN4 contacts", dx=-0.055)
axE = fig.add_subplot(gs[2, 1]); imgpanel(axE, "E_light.png", "E", "VL-CLDN4 contacts", dx=-0.055)

assert_no_text_overlap(fig, "Figure2_composite")
fig.savefig(F + "Figure2_composite.pdf", bbox_inches="tight")
fig.savefig(F + "Figure2_composite.png", dpi=600, bbox_inches="tight")
print("  wrote Figure2_composite")
