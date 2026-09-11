#!/usr/bin/env python3
"""
Figure 2B: STRIVE specificity landscape, regenerated as vector from the deposited
scoring table (replaces the raster panel from
notebooks/scfv_specificity_decision_notebook.ipynb).

Reproduces the notebook's logic exactly - R_spec = STRIVE(CLDN4) / sqrt(mean_off x max_off),
tiered at 1.0 / 1.2 / 1.5. The notebook's category labels used emoji that render as missing
glyphs; they are plain text here.

An earlier version paired this with a 39-row bar chart of the worst-case ratio R_max. That
chart has been dropped: 39 legible rows need about 135 mm of height, and the panel gets
roughly a third of that in any realistic layout, so its labels came out at 3 pt. The
worst-case result it carried is reported in Table 1 (per-candidate R_max and a pass/fail
column), in Section 3.4, in the Figure 2 legend, and as the annotation below - none of which
depends on reading 39 tick labels.

SIZE. Drawn 90 x 75 mm so that placing it at 90 mm wide reproduces the type sizes here 1:1
(7.5 pt body, 6.8 pt legend). Placed narrower, multiply through: at 70 mm the body type is
7.5 x 70/90 = 5.8 pt, below the legibility floor.

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
MM = 1 / 25.4
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
# ordered tiers -> one-hue sequential ramp (the quantity is ordinal, not nominal)
TIERS = [("Not specific (<= 1.0)", "#b7d3f6"), ("Acceptable (1.0-1.2)", "#6da7ec"),
         ("Good (1.2-1.5)", "#2a78d6"), ("Excellent (> 1.5)", "#104281")]
plt.rcParams.update({"font.size": 7.5, "axes.edgecolor": "#c9c8c4", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF})

d = pd.read_csv(O + "TableS_specificity_full.csv")
d["label"] = d.scFv.str.replace("scfv_", "scFv-", regex=False).str.replace("_", "-", regex=False)
d["tier"] = d.R_spec.map(lambda r: 3 if r > 1.5 else 2 if r > 1.2 else 1 if r > 1.0 else 0)
d["colour"] = d.tier.map(lambda t: TIERS[t][1])

fig, ax = plt.subplots(figsize=(90 * MM, 75 * MM))
ax.axhline(1.0, color="#9a9994", lw=0.8, ls="--")
ax.scatter(d.CLDN4, d.R_spec, s=30, c=d.colour, edgecolor=SURF, linewidth=0.6, zorder=3)
# scFv-48-1 sits just right of and below scFv-73-0, so the default up-left offset
# would drop its label on top of the scFv-73-0 marker; it is placed below instead.
OFFSET = {"scFv-48-1": (-6, -11)}
for _, r in d[d.R_spec > 1.2].iterrows():
    ax.annotate(r.label, (r.CLDN4, r.R_spec), textcoords="offset points",
                xytext=OFFSET.get(r.label, (-6, 6)), ha="right", fontsize=6.8, color=INK)
ax.set_xlabel("On-target CLDN4 STRIVE score")
ax.set_ylabel("R$_{spec}$  (specificity ratio)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor=c,
                          markeredgecolor=SURF, label=l) for l, c in TIERS],
          fontsize=6.8, frameon=False, loc="upper left", handletextpad=0.4,
          borderpad=0.2, labelspacing=0.35)
# the worst-case result the dropped bar chart carried, stated rather than plotted
ax.text(0.985, 0.03, "No candidate clears the prespecified\nworst-case margin "
        "(R$_{max}$ $\\geq$ 1.20; Table 1)", transform=ax.transAxes, fontsize=6.8,
        color="#eb6834", ha="right", va="bottom", linespacing=1.35)
for s in ("top", "right"): ax.spines[s].set_visible(False)

fig.tight_layout()
assert_no_text_overlap(fig, "Figure2B_specificity")
fig.savefig(F + "Figure2B_specificity.pdf", bbox_inches="tight")
fig.savefig(F + "Figure2B_specificity.png", dpi=600, bbox_inches="tight")
print(f"  {len(d)} candidates; R_spec > 1.5: {(d.R_spec>1.5).sum()}; "
      f"R_max >= 1.2: {(d.R_max>=1.2).sum()}")
