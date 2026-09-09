#!/usr/bin/env python3
"""
Supplementary Figure 3 -- STRIVE scores for clinical-stage claudin-targeting
monoclonal antibodies across the 25-antigen claudin panel.

Regenerated for the revision to:
  (i)  collapse eclutatug and lixudebart into a single row. These are two distinct
       Alentis anti-CLDN1 antibodies (ALE.C04 and ALE.F02) that share a BYTE-IDENTICAL
       Fv - identical VH, identical VL, identical CDRs - and differ only by the IgG1
       Fc substitutions L234F/L235E/P331S, which lie outside the modelled
       antigen-binding interface. Their identical AlphaFold3 metrics are therefore the
       correct deterministic result, not a data-entry error, and the pair doubles as an
       internal control on pipeline determinism. The WHO INN stems encode exactly this
       relationship: "-tug" denotes an unmodified immunoglobulin and "-bart" one with
       engineered constant regions.
  (ii) replace the original red-yellow-green diverging colour map, which encoded a pure
       magnitude on a diverging scale, with a single-hue sequential ramp.
Cognate targets are taken from the pipeline source
(DesCARTes_ozan/scripts/mab_offtarget_af3_combinator.py, MAB_ON_TARGETS) and each
antibody's cognate cell is outlined so that identity is never carried by colour alone.

Usage:  python3 R8_supp_fig3_mab_heatmap.py
Outputs: figures/SupplementaryFigure3_mAb_panel.pdf  (vector, for submission)
         figures/SupplementaryFigure3_mAb_panel.png  (600 dpi preview)
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

SRC = "/mnt/ssd1/Projects/DesCARTes_wd/DesCARTes_ozan/results/mab_offtarget_scores.csv"
OUTDIR = "figures/"
# Fv-identical pair collapsed into one row (see module docstring)
MERGE = ("eclutatug", "lixudebart")
MERGE_LABEL = "eclutatug / lixudebart"
# "sonesitamab" in the source tables is a misspelling of the INN sonesitatug
# (sonesitatug vedotin, CMG901/AZD0901; KYM Biosciences/AstraZeneca), an
# anti-CLDN18.2 ADC. Relabelled here; the source key is unchanged.
RENAME = {"sonesitamab": "sonesitatug vedotin"}

# Cognate targets, from the pipeline source (MAB_ON_TARGETS)
ON_TARGET = {
    "sonesitatug vedotin": "CLDN18.2",
    MERGE_LABEL: "CLDN1", "lixudebart": "CLDN1", "eclutatug": "CLDN1", "ixotatug": "CLDN6",
    "garetatug": "CLDN18.2", "omectatug": "CLDN18.2", "osemitamab": "CLDN18.2",
    "sonesitamab": "CLDN18.2", "tecotabart": "CLDN18.2",
    "gresonitamab": "CLDN18.2", "zolbetuximab": "CLDN18.2",
}

# Sequential ramp: one hue, light -> dark (blue 100 -> 700). A magnitude in [0,1]
# takes a sequential ramp, never a diverging or rainbow map.
BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
        "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue", BLUE)
SURFACE, INK, INK_2 = "#fcfcfb", "#0b0b0b", "#52514e"

def col_label(c):
    c = c.replace("cld", "CLDN").replace("_", ".")
    return c

def col_key(c):
    b = c.replace("CLDN", "")
    return (int(b.split(".")[0]), int(b.split(".")[1]) if "." in b else 0)

m = pd.read_csv(SRC)
piv = m.pivot_table(index="MAb", columns="Target", values="Scientific_Score")

# verify the pair really is identical before collapsing, then collapse
a, b = MERGE
assert a in piv.index and b in piv.index, "merge pair not present"
assert np.allclose(piv.loc[a].values, piv.loc[b].values, atol=1e-9), \
    "MERGE pair is no longer identical - re-check before collapsing"
piv = piv.drop(index=[b]).rename(index={a: MERGE_LABEL}).rename(index=RENAME)
print(f"collapsed {a} + {b} -> '{MERGE_LABEL}' (verified identical across all antigens)")
piv.columns = [col_label(c) for c in piv.columns]
piv = piv[sorted(piv.columns, key=col_key)]

# rows grouped by cognate target, then by cognate-target score (best first)
piv["_t"] = [ON_TARGET.get(i, "?") for i in piv.index]
piv["_s"] = [piv.loc[i, ON_TARGET.get(i)] if ON_TARGET.get(i) in piv.columns else -1
             for i in piv.index]
piv = piv.sort_values(["_t", "_s"], ascending=[True, False])
targets = piv.pop("_t").to_dict()
piv.pop("_s")

nrow, ncol = piv.shape
fig, ax = plt.subplots(figsize=(0.46 * ncol + 4.4, 0.46 * nrow + 2.6))
fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
norm = Normalize(0.0, 1.0)
im = ax.imshow(piv.values, cmap=CMAP, norm=norm, aspect="equal")

# 2px surface gap between cells
ax.set_xticks(np.arange(-0.5, ncol, 1), minor=True)
ax.set_yticks(np.arange(-0.5, nrow, 1), minor=True)
ax.grid(which="minor", color=SURFACE, linewidth=2)
ax.tick_params(which="minor", length=0)

for i in range(nrow):
    for j in range(ncol):
        v = piv.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.0,
                color="#ffffff" if v > 0.62 else INK)

# outline each antibody's cognate cell -- identity is not colour-alone
cols = list(piv.columns)
for i, ab in enumerate(piv.index):
    t = targets.get(ab)
    if t in cols:
        ax.add_patch(Rectangle((cols.index(t) - 0.5, i - 0.5), 1, 1, fill=False,
                               edgecolor="#eb6834", linewidth=2.0, zorder=5))

ax.set_xticks(range(ncol)); ax.set_yticks(range(nrow))
ax.set_xticklabels(cols, rotation=90, fontsize=7.5, color=INK)
ax.set_yticklabels([f"{ab}  ({targets.get(ab,'?')})" for ab in piv.index],
                   fontsize=8, color=INK)
ax.set_xlabel("Claudin family antigen", fontsize=9, color=INK_2, labelpad=10)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)

cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.015)
cb.set_label("STRIVE score (model-derived prioritization value)", fontsize=8, color=INK_2)
cb.ax.tick_params(labelsize=7, length=2, color=INK_2)
cb.outline.set_visible(False)

ax.set_title("Clinical-stage claudin-targeting monoclonal antibodies scored against "
             "the 25-antigen claudin panel",
             fontsize=10, color=INK, pad=14, loc="left")
fig.text(0.005, 0.012,
         "Orange outline marks each antibody's cognate target. Scores are model-derived and are "
         "not measured affinities; differences below ~0.07 units lie within the run-to-run\n"
         "variability of the predictor. eclutatug (ALE.C04) and lixudebart (ALE.F02) share an "
         "identical Fv and differ only in the IgG1 Fc, outside the modelled interface, so they\n"
         "score identically and are shown as one row.",
         fontsize=7, color=INK_2, va="bottom")
fig.tight_layout(rect=[0, 0.105, 1, 0.94])
fig.savefig(OUTDIR + "SupplementaryFigure3_mAb_panel.pdf", facecolor=SURFACE)
fig.savefig(OUTDIR + "SupplementaryFigure3_mAb_panel.png", dpi=600, facecolor=SURFACE)
print(f"wrote figures for {nrow} antibodies x {ncol} antigens")

# ---- how often is the cognate target the top-ranked antigen? -----------------
rows = []
for ab in piv.index:
    t = targets.get(ab)
    s = piv.loc[ab]
    rank = int((s > s[t]).sum()) + 1 if t in cols else None
    rows.append(dict(mAb=ab, cognate_target=t,
                     cognate_score=round(s[t], 3) if t in cols else None,
                     top_antigen=s.idxmax(), top_score=round(s.max(), 3),
                     cognate_rank_of_25=rank))
r = pd.DataFrame(rows)
r.to_csv("analysis/output/mab_panel_cognate_ranks.csv", index=False)
print(r.to_string(index=False))
print(f"\ncognate ranked first: {(r.cognate_rank_of_25 == 1).sum()} / {len(r)}")
