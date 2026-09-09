#!/usr/bin/env python3
"""
Revision analysis 5b -- LuCA CLDN4 stage comparison, reproducing the published
unit of analysis and then correcting it.
Addresses Reviewer 1 major point 2 and Reviewer 2 major point 4.

The submitted analysis collapsed the five normal lung epithelial cell types into
a single "Epithelial (normal lung)" pseudobulk per sample (195 profiles) and used
one malignant-cell pseudobulk per tumour sample (248 profiles), giving 443
profiles, and compared groups with an unpaired Mann-Whitney U test without
multiplicity correction.

Here we (a) reproduce that design exactly, (b) add Benjamini-Hochberg FDR
correction over the five stage comparisons, (c) repeat the analysis with one
profile per donor to remove pseudoreplication across samples from the same
individual, and (d) test the monotonic-progression claim explicitly.

Input : output/luca_cldn4_profiles.csv (per sample x cell_type raw sums,
        produced by R5_luca_pseudobulk_stats.py)
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

STAGES = ["I", "II", "III", "III or IV", "IV"]
REF = "Epithelial (normal lung)"
c = pd.read_csv("output/luca_cldn4_profiles.csv")

def agg(df, keys):
    g = df.groupby(keys, observed=True)[["g", "tot"]].sum().reset_index()
    g["expr"] = np.log1p(g.g / g.tot * 1e6)
    return g

# (A) published design: one profile per sample per GROUP (epithelial types collapsed)
sample_lvl = agg(c, ["sample", "donor_id", "study", "group"])
# (B) donor-level: one profile per donor per group
donor_lvl = agg(c, ["donor_id", "study", "group"])

lines = []
def emit(s=""):
    print(s); lines.append(s)

emit("=" * 92)
emit("LuCA PSEUDOBULK CLDN4 -- REPRODUCTION AND CORRECTION")
emit("=" * 92)
for lab, t in [("sample-level (published design)", sample_lvl), ("donor-level", donor_lvl)]:
    n = t.group.value_counts().reindex([REF] + STAGES).fillna(0).astype(int)
    emit(f"\n{lab}: {len(t)} profiles")
    emit("  " + "; ".join(f"{k} n={v}" for k, v in n.items()))

rows = []
for lab, t in [("sample-level (published design)", sample_lvl), ("donor-level", donor_lvl)]:
    r = t.loc[t.group == REF, "expr"].values
    for s in STAGES:
        v = t.loc[t.group == s, "expr"].values
        if not len(v):
            continue
        u = stats.mannwhitneyu(v, r, alternative="two-sided")
        # rank-biserial effect size
        rb = 2 * u.statistic / (len(v) * len(r)) - 1
        rows.append(dict(unit=lab, stage=s, n=len(v), n_ref=len(r),
                         median=np.median(v), median_ref=np.median(r),
                         log2FC=np.median(v) - np.median(r),
                         U=u.statistic, rank_biserial=rb, p_raw=u.pvalue))
res = pd.DataFrame(rows)
for lab, sub in res.groupby("unit"):
    res.loc[sub.index, "p_BH"] = multipletests(sub.p_raw, method="fdr_bh")[1]
res["signif_BH"] = np.where(res.p_BH < 0.001, "***",
                    np.where(res.p_BH < 0.01, "**",
                      np.where(res.p_BH < 0.05, "*", "ns")))
res.to_csv("output/luca_stage_stats_FINAL.csv", index=False)

emit("")
emit(res.round(4).to_string(index=False))
emit("")

for lab, t in [("sample-level", sample_lvl), ("donor-level", donor_lvl)]:
    tt = t[t.group.isin(STAGES)].copy()
    tt["rank"] = tt.group.map({s: i for i, s in enumerate(STAGES)})
    rho, pr = stats.spearmanr(tt["rank"], tt.expr)
    kw = stats.kruskal(*[tt.loc[tt.group == s, "expr"].values for s in STAGES])
    emit(f"[{lab}] monotonic trend across tumour stages: Spearman rho = {rho:.3f} "
         f"(p = {pr:.3g}); Kruskal-Wallis H = {kw.statistic:.2f} (p = {kw.pvalue:.3g})")

emit("")
emit("Median log1p(CPM) by group (sample-level):")
emit(sample_lvl.groupby("group").expr.median().reindex([REF] + STAGES).round(3).to_string())
emit("")
emit("Cohorts contributing to each group (donor-level):")
comp = pd.crosstab(donor_lvl.group, donor_lvl.study)
emit("  " + "; ".join(f"{g}: {(comp.loc[g] > 0).sum()} studies"
                      for g in [REF] + STAGES if g in comp.index))
emit("  NOTE: the 'III or IV' group is contributed entirely by a single cohort")
emit("  (Wu_Zhou_2021), so that comparison is confounded with study of origin.")

open("output/luca_report_FINAL.txt", "w").write("\n".join(lines) + "\n")
print("\nSaved -> output/luca_stage_stats_FINAL.csv, luca_report_FINAL.txt")
