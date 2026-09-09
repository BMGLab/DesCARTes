#!/usr/bin/env python3
"""
Regenerate the figure panels that can be rebuilt from deposited data, as vector PDF
plus 600 dpi PNG. Addresses the Editorial Office resolution request and Reviewer 2's
presentation points.

Panels produced:
  Fig1B  LuCA donor-level pseudobulk CLDN4 by stage, FDR-corrected
  Fig1C  HPA IHC scoring, normal lung vs NSCLC
  Fig1D  Matched cohort IHC scoring (paired)
  Fig3   MD: global RMSD, Rg, and the new interface-resolved panel
  FigS2  RoseTTAFold2 filtering metrics for all 1,000 designs
  FigS4  STRIVE cross-reactivity landscape, 39 candidates x 25 antigens
  FigS5  CLDN4 across normal human tissues (new)

Figure 1A (TCGA bulk) and Figure 2 (workflow schematic + structure renders) are not
regenerated here: the TCGA expression matrix is not in the deposit and Figure 2 is
hand-assembled artwork.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcheck import assert_no_text_overlap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

O, F = "output/", "../figures/"
BLUE = ["#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5",
        "#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue", BLUE)
S1, S2c, S3c = "#2a78d6", "#eb6834", "#1baf7a"          # categorical slots 1-3
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": "#c9c8c4",
                     "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF})

def save(fig, name):
    assert_no_text_overlap(fig, name)
    fig.savefig(F + name + ".pdf", bbox_inches="tight")
    fig.savefig(F + name + ".png", dpi=600, bbox_inches="tight")
    plt.close(fig); print("  wrote", name)

def bare(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)

# ---------------------------------------------------------------- Fig 1B ----
prof = pd.read_csv(O + "luca_cldn4_profiles.csv")
don = prof.groupby(["donor_id", "group"], observed=True)[["g", "tot"]].sum().reset_index()
don["expr"] = np.log1p(don.g / don.tot * 1e6)
ORD = ["Epithelial (normal lung)", "I", "II", "III", "III or IV", "IV"]
LAB = ["Normal\nepithelium", "I", "II", "III", "III/IV", "IV"]
st = pd.read_csv(O + "luca_stage_stats_FINAL.csv")
st = st[st.unit == "donor-level"].set_index("stage")
fig, ax = plt.subplots(figsize=(5.6, 3.5))
data = [don.loc[don.group == g, "expr"].values for g in ORD]
shades = [CMAP(v) for v in np.linspace(0.18, 0.88, len(ORD))]
bp = ax.boxplot(data, patch_artist=True, widths=0.58, showfliers=False,
                medianprops=dict(color=INK, lw=1.3),
                whiskerprops=dict(color="#9a9994"), capprops=dict(color="#9a9994"),
                boxprops=dict(lw=0))
for patch, c in zip(bp["boxes"], shades):
    patch.set_facecolor(c); patch.set_edgecolor(SURF); patch.set_linewidth(2)
rng = np.random.default_rng(0)
for i, d in enumerate(data, start=1):
    ax.scatter(rng.normal(i, 0.075, len(d)), d, s=5, color=INK, alpha=0.30,
               linewidths=0, zorder=3)
ytop = max(v.max() for v in data)
for i, g in enumerate(ORD[1:], start=2):
    if g in st.index:
        sig = st.loc[g, "signif_BH"]
        ax.text(i, ytop + 0.42, sig, ha="center", va="bottom", fontsize=9,
                color=INK if sig != "ns" else INK2)
# n goes in the tick label, not inside the axes, so it cannot land on a data point
ax.set_xticks(range(1, len(ORD) + 1))
ax.set_xticklabels([f"{l}\nn={len(d)}" for l, d in zip(LAB, data)], fontsize=7.5)
ax.set_ylabel("CLDN4  log(CPM+1)", color=INK2)
ax.set_ylim(top=ytop + 1.15)
ax.set_title("Donor-level CLDN4 in the LuCA atlas (Benjamini-Hochberg corrected)",
             fontsize=8.5, color=INK, loc="left")
bare(ax); save(fig, "Figure1B_LuCA_stage")

# ------------------------------------------------------------- Fig 1C/1D ----
LV = ["Not detected", "Low", "Medium", "High"]
cols = [CMAP(v) for v in (0.10, 0.35, 0.62, 0.88)]
hpa = {"Normal lung\n(n=4)": [3, 1, 0, 0], "NSCLC\n(n=11)": [0, 2, 9, 0]}
ihc = pd.read_csv(O + "TableS_IHC_pairs.csv")
t = ihc.tumour_intensity.value_counts()
own = {"Adjacent normal\n(n=45)": [45, 0, 0, 0],
       "NSCLC tumour\n(n=45)": [0, int(t.get(1, 0)), int(t.get(2, 0)), int(t.get(3, 0))]}
for name, dat, title, note in [
    ("Figure1C_HPA_IHC", hpa, "CLDN4 IHC, Human Protein Atlas", ""),
    ("Figure1D_matched_IHC", own, "CLDN4 IHC, matched cohort",
     "all 45 pairs concordant; exact Wilcoxon signed-rank p = 5.7 x 10$^{-14}$")]:
    fig, ax = plt.subplots(figsize=(3.5, 3.2))
    keys = list(dat); bottoms = np.zeros(len(keys))
    for li, lv in enumerate(LV):
        vals = np.array([dat[k][li] for k in keys], float)
        pct = vals / np.array([sum(dat[k]) for k in keys]) * 100
        ax.bar(keys, pct, bottom=bottoms, color=cols[li], label=lv,
               edgecolor=SURF, linewidth=2, width=0.6)
        for xi, (p_, b_) in enumerate(zip(pct, bottoms)):
            if p_ >= 8:
                ax.text(xi, b_ + p_ / 2, f"{p_:.0f}%", ha="center", va="center",
                        fontsize=7, color="#ffffff" if li >= 2 else INK)
        bottoms += pct
    ax.set_ylabel("% of cases", color=INK2); ax.set_ylim(0, 100)
    ax.set_title(title, fontsize=8.5, color=INK, loc="left")
    if note: ax.text(0, -0.22, note, transform=ax.transAxes, fontsize=6.5, color=INK2)
    ax.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(cols, LV)],
              fontsize=6.5, frameon=False, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), title="Staining", title_fontsize=6.5)
    bare(ax); save(fig, name)

# ----------------------------------------------------------------- Fig 3 ----
ts = pd.read_csv(O + "md_interface_timeseries.csv")
bsa = pd.read_csv(O + "md_bsa.csv")
ok = ts.scFv_rmsd_antigen_aligned <= 20
fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.1))
a = axes[0]
a.plot(ts.time_ns[ok], ts.complex_backbone_rmsd[ok], lw=1.0, color=S1, label="Complex backbone")
a.plot(ts.time_ns[ok], ts.CLDN4_backbone_rmsd[ok], lw=1.0, color=S3c, label="CLDN4 backbone")
a.plot(ts.time_ns[ok], ts.interface_rmsd[ok], lw=1.0, color=S2c, label="Interface")
a.axvspan(0, 100, color="#efeee9", zorder=0)
a.set_xlabel("Time (ns)"); a.set_ylabel("RMSD (Å)")
a.legend(fontsize=6.5, frameon=False, loc="upper left")
# label the shaded band along the bottom of the axes, clear of the legend
a.text(50, a.get_ylim()[0] + 0.04 * (a.get_ylim()[1] - a.get_ylim()[0]),
       "relaxation", ha="center", va="bottom", fontsize=6, color=INK2)
a.set_title("Structural stability", fontsize=8.5, color=INK, loc="left"); bare(a)
b = axes[1]
b.plot(bsa.time_ns, bsa.bsa, lw=1.0, color=S1)
b.set_xlabel("Time (ns)"); b.set_ylabel("Buried surface area (Å$^2$)")
b.set_title("Interface size", fontsize=8.5, color=INK, loc="left"); bare(b)
c = axes[2]
c.plot(ts.time_ns, ts.atom_contacts, lw=1.0, color=S1, label="Heavy-atom contacts")
c.set_xlabel("Time (ns)"); c.set_ylabel("Contacts (<4.5 Å)", color=S1)
c.tick_params(axis="y", colors=S1)
c2 = c.twinx()
c2.plot(ts.time_ns, ts.interface_hbonds, lw=1.0, color=S2c)
c2.set_ylabel("Hydrogen bonds", color=S2c); c2.tick_params(axis="y", colors=S2c)
c2.spines["top"].set_visible(False)
c.set_title("Interface persistence", fontsize=8.5, color=INK, loc="left"); bare(c)
fig.suptitle("500 ns scFv-73-0/CLDN4 membrane simulation — the predicted pose is retained",
             fontsize=9, color=INK, x=0.005, ha="left", y=1.04)
fig.tight_layout(); save(fig, "Figure3_MD_panels")

# ---------------------------------------------------------------- Fig S2 ----
sc = pd.read_csv("/mnt/ssd1/Projects/DesCARTes_wd/DesCARTes_ozan/results/3_rf2.sc", sep="\t")
fig, axes = plt.subplots(1, 4, figsize=(12.2, 2.9))
for ax, (col, lab, thr, keep_below) in zip(axes, [
        ("pae", "Global PAE", 10, True), ("pred_lddt", "Predicted lDDT", None, None),
        ("framework_aligned_cdr_rmsd", "Framework-aligned CDR RMSD (Å)", None, None),
        ("target_aligned_antibody_rmsd", "Target-aligned antibody RMSD (Å)", 5, True)]):
    ax.hist(sc[col].dropna(), bins=45, color=CMAP(0.55), edgecolor=SURF, linewidth=0.4)
    if thr is not None:
        ax.axvline(thr, color=S2c, lw=1.4, ls="--")
        lo, hi = ax.get_xlim()
        ax.axvspan(lo, thr, color=S2c, alpha=0.09)
        n = int((sc[col] <= thr).sum())
        ax.text(0.97, 0.93, f"≤{thr}: n={n}", transform=ax.transAxes, ha="right",
                fontsize=6.5, color=INK2)
    ax.set_xlabel(lab); ax.set_ylabel("Designs" if ax is axes[0] else ""); bare(ax)
fig.suptitle(f"RoseTTAFold2 filtering metrics for all {len(sc):,} de novo designs",
             fontsize=9, color=INK, x=0.005, ha="left", y=1.03)
fig.tight_layout(); save(fig, "SupplementaryFigure2_RF2_metrics")

# ---------------------------------------------------------------- Fig S4 ----
mat = pd.read_csv(O + "STRIVE_matrix_39x25.csv", index_col=0)
def ckey(c):
    b = c.replace("CLDN", ""); return (int(b.split(".")[0]), int(b.split(".")[1]) if "." in b else 0)
mat = mat[sorted(mat.columns, key=ckey)]
mat.index = [i.replace("scfv_", "scFv-").replace("_", "-") for i in mat.index]
mat = mat.loc[mat["CLDN4"].sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(9.6, 9.0))
im = ax.imshow(mat.values, cmap=CMAP, vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=6.5)
ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=6.5)
ax.add_patch(plt.Rectangle((list(mat.columns).index("CLDN4") - 0.5, -0.5), 1, mat.shape[0],
                           fill=False, edgecolor=S2c, lw=2))
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.012)
cb.set_label("STRIVE score", fontsize=7, color=INK2); cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=6.5)
ax.set_title("Predicted cross-reactivity landscape: 39 candidates × 25 antigens\n"
             "(CLDN4 on-target column outlined)", fontsize=9, color=INK, loc="left")
save(fig, "SupplementaryFigure4_crossreactivity")

# ---------------------------------------------------------------- Fig S5 ----
nt = pd.read_csv(O + "TableS_CLDN4_normal_tissue.csv")
NUM = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
nt["v"] = nt.level.map(NUM)
piv = (nt.groupby(["organ_system", "tissue"]).v.max().reset_index()
         .sort_values(["organ_system", "v"], ascending=[True, False]))
fig, ax = plt.subplots(figsize=(6.4, 9.4))
ypos = np.arange(len(piv))[::-1]
ax.barh(ypos, piv.v, color=[CMAP(0.12 + 0.26 * v) for v in piv.v],
        edgecolor=SURF, linewidth=1.6, height=0.74)
ax.set_yticks(ypos); ax.set_yticklabels(piv.tissue, fontsize=6.8)
ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(list(NUM), fontsize=7)
ax.set_xlabel("Maximum CLDN4 staining in any assessed cell type", color=INK2)
prev = None
for y, sysname in zip(ypos, piv.organ_system):
    if sysname != prev:
        ax.text(3.12, y, sysname, fontsize=6.2, color=INK2, va="center"); prev = sysname
ax.set_xlim(0, 3.05)
ax.set_title("CLDN4 across 48 normal human tissues (HPA v22.0)\n"
             "42 cell types in 25 tissues stain medium or high",
             fontsize=9, color=INK, loc="left")
bare(ax); save(fig, "SupplementaryFigure5_normal_tissue")
print("done")
