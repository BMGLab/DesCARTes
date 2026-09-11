#!/usr/bin/env python3
"""
Composite Figure 2 (panels A-E).

  A  De novo design and off-target screening workflow
  B  STRIVE specificity landscape (regenerated as vector from the deposited scoring
     table, replacing the raster panel from scfv_specificity_decision_notebook.ipynb)
  C  Predicted scFv-73-0/CLDN4 complex
  D  VH-CLDN4 interface contacts
  E  VL-CLDN4 interface contacts

Panels A and C-E are the authors' renders, taken at native resolution from the embedded
images of Gocmenetal_figures_highquality.pdf and composited onto white using their
accompanying soft masks (which is how the black PyMOL background is knocked out in the
published figure). Effective resolution in the published layout is 506 ppi for panel C
and ~320 ppi for D and E, i.e. already above the 300 dpi print threshold - the low ppi
reported by the standalone CLDN4_scfv73_0*.pdf wrappers reflects their 960x540 pt page
size, not the images.

Outputs: figures/Figure2_composite.{pdf,png}
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
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
# Two defects in the source artwork are corrected as overlays (data coordinates of the
# 2916x1944 raster) rather than by editing the authors' file: the design step was
# labelled "RFDiffussion", and the MD counter-ions were drawn as Na+ although the
# production system was neutralised with 61 K+ and 45 Cl- (~0.15 M KCl; Supplementary Methods S5).
axA.add_patch(Rectangle((812, 649), 259, 54, facecolor="#edf7f9", ec="none", zorder=4))
axA.text(941, 676, "RFdiffusion", fontsize=7.5, color=INK, ha="center", va="center", zorder=5)
for cx, cy in [(2588, 1520), (2426, 1617), (2703, 1748)]:
    axA.add_patch(Circle((cx, cy), 18.5, facecolor="#288a7a", ec="none", zorder=4))
    axA.text(cx, cy, "K", fontsize=3.4, color="white", ha="center", va="center", zorder=5)

# ------------------------------------------------------------------- B ------
axB = fig.add_subplot(gs[1, 0])
d = pd.read_csv(O + "TableS_specificity_full.csv")
d["label"] = d.scFv.str.replace("scfv_", "scFv-", regex=False).str.replace("_", "-", regex=False)
d["tier"] = d.R_spec.map(lambda r: 3 if r > 1.5 else 2 if r > 1.2 else 1 if r > 1.0 else 0)
axB.axhline(1.0, color="#9a9994", lw=0.8, ls="--")
axB.scatter(d.CLDN4, d.R_spec, s=26, c=d.tier.map(lambda t: TIERS[t][1]),
            edgecolor=SURF, linewidth=0.5, zorder=3)
OFFSET = {"scFv-48-1": (-5, -10)}     # keeps its label off the scFv-73-0 marker
for _, r in d[d.R_spec > 1.2].iterrows():
    axB.annotate(r.label, (r.CLDN4, r.R_spec), textcoords="offset points",
                 xytext=OFFSET.get(r.label, (-5, 5)), ha="right", fontsize=6.8, color=INK)
axB.set_xlabel("On-target CLDN4 STRIVE score")
axB.set_ylabel("R$_{spec}$  (specificity ratio)")
axB.legend(handles=[Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor=c,
                           markeredgecolor=SURF, label=l) for l, c in TIERS],
           fontsize=6.8, frameon=False, loc="upper left", handletextpad=0.4,
           borderpad=0.15, labelspacing=0.3)
axB.text(0.985, 0.03, "No candidate clears the prespecified\nworst-case margin "
         "(R$_{max}$ $\\geq$ 1.20; Table 1)", transform=axB.transAxes, fontsize=6.3,
         color="#eb6834", ha="right", va="bottom", linespacing=1.35)
for s in ("top", "right"): axB.spines[s].set_visible(False)
panel(axB, "B", dx=-0.16)

# ------------------------------------------------------------------- C ------
axC = fig.add_subplot(gs[1, 1]); imgpanel(axC, "C_complex.png", "C", dx=-0.02)
# The chain labels are vector text on the source page, not part of the embedded raster,
# so they are re-set here. Positions are the source page's own label boxes mapped into
# the padded image (the raster was padded 14% left / 20% right to give them the same
# margins they had on the page); given in data coordinates so they track the image
# rather than the letterboxed axes box.
cW, cH = mpimg.imread(F + "source_figure2/C_complex.png").shape[1::-1]
C_LABELS = [(0.033, 0.784, "Light chain", "left"),
            (0.987, 0.690, "Heavy chain", "right"),
            (0.714, 0.523, "scFv-73-0", "left"),
            (0.039, 0.251, "CLDN4", "left")]
for fx, fy, txt, ha in C_LABELS:
    axC.text(fx * cW, (1 - fy) * cH, txt, fontsize=7.5, color=INK, ha=ha, va="center")

# ---------------------------------------------------------------- D and E ---
axD = fig.add_subplot(gs[2, 0]); imgpanel(axD, "D_heavy.png", "D", "VH-CLDN4 contacts", dx=-0.055)
axE = fig.add_subplot(gs[2, 1]); imgpanel(axE, "E_light.png", "E", "VL-CLDN4 contacts", dx=-0.055)

# ------------------------------------------------- individual panel files ---
# Emitted so the figure can also be hand-assembled: vector PDF for the plotted panel,
# PNG for the renders. A and C are re-rendered rather than copied, so the stand-alone
# files carry the same corrections and labels the composite does - copying the raw
# source would silently drop them. Type is sized in source pixels, so it keeps its
# proportion to the baked-in artwork however the panel is later scaled.
import shutil
P = F + "panels/"
os.makedirs(P, exist_ok=True)
DPI = 300

def px_fontsize(cap_px, dpi=DPI):
    """Point size whose cap height is cap_px pixels of the source raster."""
    return cap_px / dpi * 72 / 0.72

def standalone(path, out, overlay):
    img = mpimg.imread(F + "source_figure2/" + path)
    h, w = img.shape[:2]
    f = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = f.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.set_axis_off()
    overlay(ax, w, h)
    f.savefig(P + out, dpi=DPI, facecolor="white")
    plt.close(f)
    # flatten to RGB: matplotlib writes RGBA, and a fully opaque alpha channel still
    # becomes an /SMask in the placed PDF, which prepress checks flag for no reason
    im = Image.open(P + out)
    if im.mode != "RGB":
        Image.alpha_composite(Image.new("RGBA", im.size, "white"), im).convert("RGB").save(P + out)

def _a_overlay(ax, w, h):
    ax.add_patch(Rectangle((812, 649), 259, 54, facecolor="#edf7f9", ec="none", zorder=4))
    ax.text(941, 676, "RFdiffusion", fontsize=px_fontsize(31), color=INK,
            ha="center", va="center", zorder=5)
    for cx, cy in [(2588, 1520), (2426, 1617), (2703, 1748)]:
        ax.add_patch(Circle((cx, cy), 18.5, facecolor="#288a7a", ec="none", zorder=4))
        ax.text(cx, cy, "K", fontsize=px_fontsize(14), color="white",
                ha="center", va="center", zorder=5)

# Panel C label size. 56 px of cap height is the largest that still clears the molecule
# on both sides at this crop - bigger labels collide with it or run off the panel - so the
# panel has to carry its own legibility through the width it is given in the layout:
#   label type size (pt) = 56 * (placed_width_pt / 1872) / 0.72
# which needs a placed width of >= 168 pt (59 mm) to reach 7 pt. At the 34 mm the current
# hand assembly gives panel C, the labels land at 4.2 pt.
C_LABEL_CAP_PX = 56

def _c_overlay(ax, w, h):
    for fx, fy, txt, ha in C_LABELS:
        ax.text(fx * w, (1 - fy) * h, txt, fontsize=px_fontsize(C_LABEL_CAP_PX),
                color=INK, ha=ha, va="center")

standalone("A_workflow.png", "Fig2_A_workflow.png", _a_overlay)
standalone("C_complex.png", "Fig2_C_complex.png", _c_overlay)
for src, dst in [("D_heavy.png", "Fig2_D_VH_contacts.png"),
                 ("E_light.png", "Fig2_E_VL_contacts.png")]:
    shutil.copyfile(F + "source_figure2/" + src, P + dst)
shutil.copyfile(F + "Figure2B_specificity.pdf", P + "Fig2_B_specificity.pdf")
print("  wrote panels/Fig2_A..E (A and C carry the corrections and labels)")

assert_no_text_overlap(fig, "Figure2_composite")
fig.savefig(F + "Figure2_composite.pdf", bbox_inches="tight")
fig.savefig(F + "Figure2_composite.png", dpi=600, bbox_inches="tight")
print("  wrote Figure2_composite")
