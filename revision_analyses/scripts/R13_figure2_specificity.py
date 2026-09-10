#!/usr/bin/env python3
"""
Figure 2B: STRIVE specificity landscape, regenerated as vector from the deposited
scoring table (replaces the raster panel from
notebooks/scfv_specificity_decision_notebook.ipynb).

Reproduces the notebook's logic exactly - R_spec = STRIVE(CLDN4) / sqrt(mean_off x max_off),
tiered at 1.0 / 1.2 / 1.5 - and additionally plots the worst-case ratio R_max, which is
the quantity relevant to safety and which the original panel did not show. The notebook's
category labels used emoji that render as missing glyphs; they are plain text here.

Outputs: figures/Figure2B_specificity.{pdf,png}
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figcheck import assert_no_text_overlap

O, F = "output/", "../figures/"
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
# ordered tiers -> one-hue sequential ramp (the quantity is ordinal, not nominal)
TIERS = [("Not specific (<= 1.0)", "#b7d3f6"), ("Acceptable (1.0-1.2)", "#6da7ec"),
         ("Good (1.2-1.5)", "#2a78d6"), ("Excellent (> 1.5)", "#104281")]
plt.rcParams.update({"font.size":8,"axes.edgecolor":"#c9c8c4","axes.labelcolor":INK2,
                     "xtick.color":INK2,"ytick.color":INK2,"figure.facecolor":SURF,
                     "axes.facecolor":SURF,"savefig.facecolor":SURF})

d = pd.read_csv(O + "TableS_specificity_full.csv")
d["label"] = d.scFv.str.replace("scfv_", "scFv-", regex=False).str.replace("_", "-", regex=False)
def tier(r):
    return 3 if r > 1.5 else 2 if r > 1.2 else 1 if r > 1.0 else 0
d["tier"] = d.R_spec.map(tier)
d["colour"] = d.tier.map(lambda t: TIERS[t][1])

fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 4.0),
                               gridspec_kw={"width_ratios": [1.0, 1.15]})

# ---- left: on-target score vs specificity ratio -----------------------------
axL.axhline(1.0, color="#9a9994", lw=0.8, ls="--")
axL.scatter(d.CLDN4, d.R_spec, s=34, c=d.colour, edgecolor=SURF, linewidth=0.6, zorder=3)
for _, r in d[d.R_spec > 1.2].iterrows():
    axL.annotate(r.label, (r.CLDN4, r.R_spec), textcoords="offset points",
                 xytext=(-6, 6), ha="right", fontsize=6.5, color=INK)
axL.set_xlabel("On-target CLDN4 STRIVE score")
axL.set_ylabel("R$_{spec}$  (specificity ratio)")
axL.set_title("Specificity ratio vs on-target score", fontsize=8.5, color=INK, loc="left")
axL.legend(handles=[Line2D([], [], marker="o", ls="", markersize=5.5, markerfacecolor=c,
                           markeredgecolor=SURF, label=l) for l, c in TIERS],
           fontsize=6.3, frameon=False, loc="upper left", handletextpad=0.4,
           borderpad=0.2, labelspacing=0.35)
for s in ("top", "right"): axL.spines[s].set_visible(False)

# ---- right: worst-case margin, which the original panel omitted -------------
s = d.sort_values("R_max")
y = np.arange(len(s))
axR.barh(y, s.R_max, color=s.colour, edgecolor=SURF, linewidth=0.5, height=0.74)
axR.axvline(1.0, color="#9a9994", lw=0.9, ls="--")
axR.axvline(1.2, color="#eb6834", lw=1.0, ls="--")
axR.text(1.205, len(s) - 0.4, "prespecified\n1.20-fold margin", fontsize=6.2,
         color="#eb6834", va="top", ha="left")
axR.set_yticks(y); axR.set_yticklabels(s.label, fontsize=5.6)
axR.set_xlabel("R$_{max}$  (on-target / worst off-target)")
axR.set_title("Worst-case margin: no candidate clears 1.20", fontsize=8.5, color=INK, loc="left")
axR.set_ylim(-0.8, len(s) - 0.2)
for sp in ("top", "right"): axR.spines[sp].set_visible(False)

fig.tight_layout()
assert_no_text_overlap(fig, "Figure2B_specificity")
fig.savefig(F + "Figure2B_specificity.pdf", bbox_inches="tight")
fig.savefig(F + "Figure2B_specificity.png", dpi=600, bbox_inches="tight")
print(f"  {len(d)} candidates; R_spec > 1.5: {(d.R_spec>1.5).sum()}; "
      f"R_max >= 1.2: {(d.R_max>=1.2).sum()}")
