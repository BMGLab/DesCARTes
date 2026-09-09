#!/usr/bin/env python3
"""
Revision analysis 4 -- Paired analysis of the matched NSCLC IHC cohort.
Addresses Reviewer 1 major point 2 and Reviewer 2 major point 4.

The 45 tumour / adjacent non-tumour pairs come from the SAME patients, so the
Wilcoxon rank-sum (Mann-Whitney U) test reported in the original submission is
inappropriate.  Here the paired Wilcoxon signed-rank test and the exact sign
test are used instead, together with the raw intensity components.

Source: Table 2 of the submitted manuscript (Appendix B).
Outputs: output/ihc_paired_stats.txt, output/TableS_IHC_pairs.csv
"""
import numpy as np
import pandas as pd
from scipy import stats

# Tumour intensity scores transcribed from Table 2 (1+ = low, 2+ = medium, 3+ = high).
tumour = [2, 3, 1, 3, 2, 1, 1, 2, 2, 1, 2, 1, 3, 2, 3, 3, 1, 3, 2, 1, 2, 1, 3,
          2, 3, 1, 2, 2, 2, 1, 2, 3, 1, 3, 1, 2, 2, 1, 1, 1, 3, 3, 1, 3, 3]
# Cases are listed in the same order as Table 2 of the manuscript. Sequential
# identifiers are used here rather than the hospital pathology accession numbers,
# which carry no analytic information; the accession numbers are given in Table 2
# for anyone who needs to cross-reference a specific specimen.
patients = [f"P{i:02d}" for i in range(1, 46)]
assert len(tumour) == 45 == len(patients)
normal = [0] * 45          # all matched adjacent non-tumour sections scored negative

t = np.array(tumour); n = np.array(normal); diff = t - n

df = pd.DataFrame({"patient": patients, "tumour_intensity": t,
                   "adjacent_normal_intensity": n, "difference": diff})
df.to_csv("output/TableS_IHC_pairs.csv", index=False)

lines = []
def p(s=""):
    print(s); lines.append(s)

p("=" * 74)
p("MATCHED NSCLC IHC COHORT -- PAIRED ANALYSIS (n = 45 patients)")
p("=" * 74)
p(f"Tumour intensity distribution : 1+ n={int((t==1).sum())} "
  f"({(t==1).mean()*100:.0f}%), 2+ n={int((t==2).sum())} ({(t==2).mean()*100:.0f}%), "
  f"3+ n={int((t==3).sum())} ({(t==3).mean()*100:.0f}%)")
p(f"Adjacent non-tumour           : negative in {int((n==0).sum())}/45 (100%)")
p(f"Median tumour score {np.median(t):.0f} (IQR {np.percentile(t,25):.0f}-"
  f"{np.percentile(t,75):.0f}); median adjacent score {np.median(n):.0f}")
p(f"Discordant pairs (tumour > adjacent): {int((diff>0).sum())}/45; "
  f"ties {int((diff==0).sum())}; reversals {int((diff<0).sum())}")
p("")

# --- Paired Wilcoxon signed-rank (exact) --------------------------------------
w = stats.wilcoxon(t, n, alternative="two-sided", zero_method="wilcox", method="exact")
p(f"Wilcoxon signed-rank (exact, paired) : W = {w.statistic:.1f}, p = {w.pvalue:.3g}")

# --- Exact sign test (binomial) ----------------------------------------------
n_pos, n_nonzero = int((diff > 0).sum()), int((diff != 0).sum())
b = stats.binomtest(n_pos, n_nonzero, 0.5, alternative="two-sided")
p(f"Exact sign test                      : {n_pos}/{n_nonzero} positive, p = {b.pvalue:.3g}")

# --- McNemar on positive/negative dichotomy ----------------------------------
# every tumour positive, every matched normal negative -> b = 45, c = 0
mc_p = stats.binomtest(45, 45, 0.5, alternative="two-sided").pvalue
p(f"Exact McNemar (positive vs negative) : b = 45, c = 0, p = {mc_p:.3g}")

# --- Effect size --------------------------------------------------------------
z = stats.norm.isf(w.pvalue / 2)
p(f"Matched-pairs rank-biserial correlation r = 1.000 (all 45 pairs concordant)")
p(f"Approximate |z| = {z:.2f}")
p("")

# --- What the original (incorrect) unpaired test gave -------------------------
u = stats.mannwhitneyu(t, n, alternative="two-sided")
p(f"[for comparison] Unpaired Mann-Whitney U as originally reported : "
  f"U = {u.statistic:.1f}, p = {u.pvalue:.3g}")
p("")
p("Conclusion: the paired analyses give the same direction and remain highly")
p("significant, so the substantive claim is unchanged; only the test and the")
p("reported statistics are corrected.")

open("output/ihc_paired_stats.txt", "w").write("\n".join(lines) + "\n")
