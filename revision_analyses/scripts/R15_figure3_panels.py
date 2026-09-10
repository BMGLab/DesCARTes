#!/usr/bin/env python3
"""
Figure 3, as five panels, emitted both individually (for hand-assembly in Illustrator)
and as a composite.

  A  Simulated system: scFv-73-0/CLDN4 in an asymmetric mammalian plasma membrane
  B  RMSD - complex backbone, CLDN4 backbone, interface
  C  Radius of gyration
  D  Buried surface area                (new; Reviewer 1, point 5)
  E  Interface persistence - contacts and hydrogen bonds   (new; Reviewer 1, point 5)

The submitted Figure 3 had three panels (render, RMSD, R_g). B-E are regenerated as
vector from the deposited trajectory analysis; D and E are the interface-resolved
metrics Reviewer 1 asked for, which no panel previously showed. Panel A is the authors'
own render, taken at native resolution (2247 x 2054, 333 ppi in the published layout)
from the embedded image of Gocmenetal_figures_highquality.pdf page 4 and composited onto
white using its soft mask.

R_g is read from the GROMACS gyrate output rather than the interface time series, which
does not carry it. The t = 449 ns frame is dropped as the periodic-image artefact already
documented in Supplementary Methods S5; it is the only frame above 27.3 A and is retained
in the deposited time series.

Outputs: figures/Figure3_composite.{pdf,png} and figures/panels/Fig3_{A..E}.*
"""
import sys, os, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from figcheck import assert_no_text_overlap

O, F = "output/", "../figures/"
P = F + "panels/"
os.makedirs(P, exist_ok=True)
RG_XVG = O + "md_radius_of_gyration.xvg"   # gmx gyrate output, copied into the deposit
ARTEFACT_NS = 449.0

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
S1, S2c, S3c = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 7.5, "axes.edgecolor": "#c9c8c4", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF})

def bare(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)

def panel(ax, letter, dx=-0.16):
    ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            color=INK, va="bottom", ha="left")

# ------------------------------------------------------------------ data ----
ts = pd.read_csv(O + "md_interface_timeseries.csv")
bsa = pd.read_csv(O + "md_bsa.csv")
ok = ts.scFv_rmsd_antigen_aligned <= 20            # drops the same artefact frame

_rg = np.loadtxt(RG_XVG, comments=["#", "@"])
rg_t, rg_a = _rg[:, 0] / 1000.0, _rg[:, 1] * 10.0  # ps -> ns, nm -> A
keep = rg_t != ARTEFACT_NS
rg_t, rg_a = rg_t[keep], rg_a[keep]

# ---------------------------------------------------------------- panels ----
def draw_A(ax):
    img = mpimg.imread(F + "source_figure3/A_membrane_system.png")
    ax.imshow(img); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    h, w = img.shape[:2]
    # "VL" is part of the render; the remaining labels were vector text on the source
    # page and are re-set here. Positions are in image pixels so they track the raster.
    for fx, fy, txt, ha, col in [(0.66, 0.955, "scFv-73-0", "center", INK),
                                 (0.86, 0.845, "VH", "center", S1),
                                 (0.12, 0.560, "CLDN4", "left", INK),
                                 (0.50, 0.030, "Asymmetric mammalian plasma membrane",
                                  "center", INK2)]:
        ax.text(fx * w, (1 - fy) * h, txt, fontsize=7.5, color=col, ha=ha, va="center")

def draw_B(ax):
    ax.plot(ts.time_ns[ok], ts.complex_backbone_rmsd[ok], lw=0.9, color=S1, label="Complex backbone")
    ax.plot(ts.time_ns[ok], ts.CLDN4_backbone_rmsd[ok], lw=0.9, color=S3c, label="CLDN4 backbone")
    ax.plot(ts.time_ns[ok], ts.interface_rmsd[ok], lw=0.9, color=S2c, label="Interface")
    ax.axvspan(0, 100, color="#efeee9", zorder=0)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (Å)")
    ax.legend(fontsize=6.2, frameon=False, loc="upper left", handlelength=1.4)
    lo, hi = ax.get_ylim()
    ax.text(50, lo + 0.03 * (hi - lo), "relaxation", ha="center", va="bottom",
            fontsize=6, color=INK2)
    bare(ax)

def draw_C(ax):
    ax.plot(rg_t, rg_a, lw=0.9, color=S1)
    ax.axvspan(0, 100, color="#efeee9", zorder=0)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Radius of gyration (Å)")
    bare(ax)

def draw_D(ax):
    ax.plot(bsa.time_ns, bsa.bsa, lw=0.9, color=S1)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Buried surface area (Å$^2$)")
    bare(ax)

def draw_E(ax):
    ax.plot(ts.time_ns, ts.atom_contacts, lw=0.9, color=S1)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Heavy-atom contacts (< 4.5 Å)", color=S1)
    ax.tick_params(axis="y", colors=S1)
    ax2 = ax.twinx()
    ax2.plot(ts.time_ns, ts.interface_hbonds, lw=0.9, color=S2c)
    ax2.set_ylabel("Hydrogen bonds", color=S2c); ax2.tick_params(axis="y", colors=S2c)
    ax2.spines["top"].set_visible(False)
    bare(ax)

# ------------------------------------------------------ individual panels ----
SIZES = {"B": (3.5, 2.5), "C": (3.5, 2.5), "D": (3.5, 2.5), "E": (3.5, 2.5)}
NAMES = {"B": "rmsd", "C": "radius_of_gyration", "D": "buried_surface_area",
         "E": "interface_persistence"}
for letter, draw in [("B", draw_B), ("C", draw_C), ("D", draw_D), ("E", draw_E)]:
    fig, ax = plt.subplots(figsize=SIZES[letter])
    draw(ax)
    fig.tight_layout()
    assert_no_text_overlap(fig, f"Fig3_{letter}")
    fig.savefig(f"{P}Fig3_{letter}_{NAMES[letter]}.pdf", bbox_inches="tight")
    plt.close(fig)
shutil.copyfile(F + "source_figure3/A_membrane_system.png",
                P + "Fig3_A_membrane_system.png")
print("  wrote panels/Fig3_A (PNG) and Fig3_B..E (vector PDF)")

# -------------------------------------------------------------- composite ----
fig = plt.figure(figsize=(7.2, 7.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.95], hspace=0.42, wspace=0.34)
axA = fig.add_subplot(gs[0:2, 0]); draw_A(axA); panel(axA, "A", dx=-0.02)
axB = fig.add_subplot(gs[0, 1]);   draw_B(axB); panel(axB, "B")
axC = fig.add_subplot(gs[1, 1]);   draw_C(axC); panel(axC, "C")
axD = fig.add_subplot(gs[2, 0]);   draw_D(axD); panel(axD, "D")
axE = fig.add_subplot(gs[2, 1]);   draw_E(axE); panel(axE, "E")

assert_no_text_overlap(fig, "Figure3_composite")
fig.savefig(F + "Figure3_composite.pdf", bbox_inches="tight")
fig.savefig(F + "Figure3_composite.png", dpi=600, bbox_inches="tight")
print("  wrote Figure3_composite")

print(f"  R_g over {len(rg_t)} frames: mean {rg_a.mean():.2f} A, range "
      f"{rg_a.min():.2f}-{rg_a.max():.2f} A; 100-500 ns median "
      f"{np.median(rg_a[rg_t >= 100]):.2f} A "
      f"(IQR {np.percentile(rg_a[rg_t >= 100], 25):.2f}-"
      f"{np.percentile(rg_a[rg_t >= 100], 75):.2f})")
