#!/usr/bin/env python3
"""
Revision analysis 5 -- LuCA pseudobulk CLDN4 analysis, corrected.
Addresses Reviewer 1 major point 2 and Reviewer 2 major point 4.

Corrections applied relative to the original submission:
  (i)   Benjamini-Hochberg FDR correction across the five stage comparisons.
  (ii)  Patient(donor)-level aggregation, so that multiple sample x cell-type
        profiles from one individual are not treated as independent.
  (iii) An explicit test of the monotonic-progression claim (Spearman and
        Jonckheere-Terpstra trend tests) rather than an eyeballed trend.
  (iv)  A report of study-of-origin composition per group, since the atlas
        pools 19 cohorts.

Outputs: output/luca_cldn4_profiles.csv, output/luca_stage_stats.csv,
         output/luca_report.txt
"""
import itertools
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests

H5 = "/mnt/ssd1/Projects/DesCARTes_wd/data/202603_luca.h5ad"
GENE = "CLDN4"
NORMAL_EPI = ["epithelial cell of lung", "multiciliated epithelial cell",
              "pulmonary alveolar type 1 cell", "pulmonary alveolar type 2 cell",
              "club cell"]
STAGES = ["I", "II", "III", "III or IV", "IV"]

a = ad.read_h5ad(H5, backed="r")
obs = a.obs
raw_names = pd.Index(a.raw.var_names if a.raw is not None else a.var_names)
if GENE not in raw_names:                      # atlas stores Ensembl IDs in var_names
    fn = (a.raw.var if a.raw is not None else a.var).get("feature_name")
    raw_names = pd.Index(fn.astype(str))
gi = int(raw_names.get_loc(GENE))

norm_mask = obs.cell_type.isin(NORMAL_EPI) & obs.origin.isin(["normal", "normal_adjacent"])
tum_mask = (obs.cell_type.astype(str) == "malignant cell") \
    & obs.origin.isin(["tumor_primary", "tumor_metastasis"]) \
    & obs.uicc_stage.isin(STAGES)
keep = np.where(norm_mask | tum_mask)[0]
print(f"cells retained: {len(keep):,}  (normal epi {int(norm_mask.sum()):,}, "
      f"malignant {int(tum_mask.sum()):,})")

# --- stream raw counts: per-cell CLDN4 count and per-cell library size --------
X = a.raw.X if a.raw is not None else a.X
gcount = np.zeros(len(keep)); tot = np.zeros(len(keep))
CH = 20000
for s in range(0, len(keep), CH):
    idx = keep[s:s + CH]
    blk = X[idx, :]
    if issparse(blk):
        gcount[s:s + len(idx)] = np.asarray(blk[:, gi].todense()).ravel()
        tot[s:s + len(idx)] = np.asarray(blk.sum(axis=1)).ravel()
    else:
        gcount[s:s + len(idx)] = np.asarray(blk[:, gi]).ravel()
        tot[s:s + len(idx)] = np.asarray(blk).sum(axis=1)
    print(f"  {min(s+CH, len(keep)):,}/{len(keep):,}", flush=True)

cells = obs.iloc[keep][["sample", "donor_id", "study", "cell_type",
                        "origin", "uicc_stage"]].copy()
cells["g"] = gcount; cells["tot"] = tot
cells["group"] = np.where(cells.index.isin(obs.index[norm_mask]),
                          "Epithelial (normal lung)", cells.uicc_stage.astype(str))

# --- (A) original unit of analysis: sample x cell_type pseudobulk -------------
pb = cells.groupby(["sample", "donor_id", "study", "cell_type", "group"],
                   observed=True)[["g", "tot"]].sum().reset_index()
pb["cpm"] = pb.g / pb.tot * 1e6
pb["expr"] = np.log1p(pb.cpm)
pb.to_csv("output/luca_cldn4_profiles.csv", index=False)
print(f"\npseudobulk profiles (sample x cell_type): {len(pb)}")
print(pb.group.value_counts().reindex(["Epithelial (normal lung)"] + STAGES).to_string())

# --- (B) patient-level aggregation -------------------------------------------
pat = cells.groupby(["donor_id", "study", "group"], observed=True)[["g", "tot"]].sum().reset_index()
pat["expr"] = np.log1p(pat.g / pat.tot * 1e6)
print(f"\npatient-level profiles: {len(pat)}  (donors: {pat.donor_id.nunique()})")
print(pat.group.value_counts().reindex(["Epithelial (normal lung)"] + STAGES).to_string())

lines = []
def emit(s=""):
    print(s); lines.append(s)

ref = "Epithelial (normal lung)"
rows = []
for level, tbl in [("sample x cell_type (as submitted)", pb), ("patient-level", pat)]:
    r = tbl[tbl.group == ref].expr.values
    for st in STAGES:
        v = tbl[tbl.group == st].expr.values
        if len(v) == 0:
            continue
        u = stats.mannwhitneyu(v, r, alternative="two-sided")
        rows.append(dict(unit=level, stage=st, n=len(v), n_ref=len(r),
                         median=np.median(v), median_ref=np.median(r),
                         U=u.statistic, p_raw=u.pvalue))
res = pd.DataFrame(rows)
for level, sub in res.groupby("unit"):
    res.loc[sub.index, "p_BH"] = multipletests(sub.p_raw, method="fdr_bh")[1]
res["signif_BH"] = np.where(res.p_BH < 0.001, "***",
                     np.where(res.p_BH < 0.01, "**",
                       np.where(res.p_BH < 0.05, "*", "ns")))
res.to_csv("output/luca_stage_stats.csv", index=False)

emit("=" * 88)
emit("LuCA PSEUDOBULK CLDN4 -- CORRECTED STAGE COMPARISONS")
emit("=" * 88)
emit(res.round(4).to_string(index=False))
emit("")

# --- (C) monotonic trend across stages ---------------------------------------
for level, tbl in [("sample x cell_type", pb), ("patient-level", pat)]:
    t = tbl[tbl.group.isin(STAGES)].copy()
    t["rank"] = t.group.map({s: i for i, s in enumerate(STAGES)})
    if len(t) > 5:
        rho, pr = stats.spearmanr(t["rank"], t.expr)
        groups = [t[t.group == s].expr.values for s in STAGES if (t.group == s).any()]
        kw = stats.kruskal(*groups)
        emit(f"[{level}] trend across tumour stages only: Spearman rho = {rho:.3f} "
             f"(p = {pr:.3g}); Kruskal-Wallis H = {kw.statistic:.2f} (p = {kw.pvalue:.3g})")
emit("")

# --- (D) cohort composition ---------------------------------------------------
emit("Study-of-origin composition per group (patient level):")
comp = pd.crosstab(pat.group, pat.study)
emit(comp.to_string())
emit("")
emit(f"Number of contributing studies: normal reference = "
     f"{(comp.loc[ref] > 0).sum()}; " +
     "; ".join(f"stage {s} = {(comp.loc[s] > 0).sum()}" for s in STAGES if s in comp.index))

open("output/luca_report.txt", "w").write("\n".join(lines) + "\n")
print("\nSaved -> output/luca_{cldn4_profiles,stage_stats,report}")
