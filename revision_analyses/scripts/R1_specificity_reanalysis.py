#!/usr/bin/env python3
"""
Revision analysis 1 -- Worst-case specificity margins for de novo scFv candidates.
Addresses Reviewer 1 major point 4 and Reviewer 2 major point 2.

Inputs  : DesCARTes_ozan/results/OTR_ontarget.csv        (39 scFv vs CLDN4)
          DesCARTes_ozan/results/OTR_offtarget_test1.csv (39 scFv x 24 homologues)
Outputs : output/TableS_specificity_full.csv   -- all 39 candidates, every metric
          output/Table1_revised.csv            -- revised main-text Table 1
          output/worstcase_margins.csv         -- per-candidate margin vs every off-target
"""
import re
import numpy as np
import pandas as pd

RES = "/mnt/ssd1/Projects/DesCARTes_wd/DesCARTes_ozan/results/"
OUT = "output/"

def scfv_id(s):
    m = re.search(r"(scfv_\d+_\d+)", str(s))
    return m.group(1) if m else None

def off_target(s):
    m = re.search(r"vs_(cld\d+(?:_\d+)?)", str(s))
    return m.group(1).upper().replace("CLD", "CLDN").replace("_", ".") if m else None

on = pd.read_csv(RES + "OTR_ontarget.csv")
off = pd.read_csv(RES + "OTR_offtarget_test1.csv")
on["scFv"] = on.Sample_ID.map(scfv_id)
off["scFv"] = off.Sample_ID.map(scfv_id)
off["target"] = off.Sample_ID.map(off_target)

# Wide STRIVE matrix: rows = 39 scFv, cols = CLDN4 + 24 homologues
mat = off.pivot_table(index="scFv", columns="target", values="Scientific_Score")
mat.insert(0, "CLDN4", on.set_index("scFv").Scientific_Score)
mat.to_csv(OUT + "STRIVE_matrix_39x25.csv")

off_cols = [c for c in mat.columns if c != "CLDN4"]

rows = []
for s in mat.index:
    tgt = mat.loc[s, "CLDN4"]
    o = mat.loc[s, off_cols].dropna()
    avg, mx = o.mean(), o.max()
    # per-off-target margin: on-target score divided by that homologue's score
    margins = tgt / o
    rows.append(dict(
        scFv=s,
        CLDN4=tgt,
        AvgOff=avg,
        MaxOff=mx,
        MaxOffTarget=o.idxmax(),
        R_avg=tgt / avg,                       # = STRIVE_CLDN4 / mean(STRIVE_off)
        R_max=tgt / mx,                        # = STRIVE_CLDN4 / max(STRIVE_off)
        R_spec=tgt / np.sqrt(avg * mx),        # geometric-mean denominator
        n_off_above_0p70=int((o >= 0.70).sum()),
        n_off_within_10pct=int((o >= 0.90 * tgt).sum()),   # near-ties with on-target
        n_off_margin_lt_1p20=int((margins < 1.20).sum()),  # fails prespecified margin
        CLDN18_2=mat.loc[s, "CLDN18.2"] if "CLDN18.2" in mat.columns else np.nan,
    ))

res = pd.DataFrame(rows).sort_values("R_spec", ascending=False)

# --- Prespecified decision rule (revised: worst-case, not geometric mean) --------
# A candidate passes only if EVERY homologue is beaten by >= 1.2-fold (R_max >= 1.2)
# AND no homologue reaches an absolute STRIVE score of 0.70.
res["pass_worstcase"] = (res.R_max >= 1.20) & (res.MaxOff < 0.70)

def category(r):
    """Revised, internally consistent tiers based on the WORST-CASE ratio R_max."""
    if r > 1.50:
        return "Excellent (R_max > 1.50)"
    if r > 1.20:
        return "Good (1.20 < R_max <= 1.50)"
    if r > 1.00:
        return "Acceptable (1.00 < R_max <= 1.20)"
    return "Not specific (R_max <= 1.00)"

res["Category_Rmax"] = res.R_max.map(category)

def category_geo(r):
    """Original geometric-mean tiers, restated with unambiguous inequalities."""
    if r > 1.50:
        return "Excellent (R_spec > 1.50)"
    if r > 1.20:
        return "Good (1.20 < R_spec <= 1.50)"
    if r > 1.00:
        return "Acceptable (1.00 < R_spec <= 1.20)"
    return "Not specific (R_spec <= 1.00)"

res["Category_Rspec"] = res.R_spec.map(category_geo)

res.round(3).to_csv(OUT + "TableS_specificity_full.csv", index=False)
res.head(15).round(3).to_csv(OUT + "Table1_revised.csv", index=False)

# --- Long-format worst-case margin table for every candidate/homologue pair -----
long = mat.reset_index().melt(id_vars=["scFv", "CLDN4"], value_vars=off_cols,
                              var_name="homologue", value_name="STRIVE_off")
long["margin"] = long.CLDN4 / long.STRIVE_off
long["fails_1p20_margin"] = long.margin < 1.20
long.sort_values(["scFv", "margin"]).round(3).to_csv(OUT + "worstcase_margins.csv", index=False)

# --- Console report -------------------------------------------------------------
pd.set_option("display.width", 220)
print("=" * 78)
print("REVISED SPECIFICITY ANALYSIS  (n = 39 candidates x 25 antigens)")
print("=" * 78)
print(res.head(15)[["scFv", "CLDN4", "AvgOff", "MaxOff", "MaxOffTarget",
                    "R_avg", "R_max", "R_spec", "Category_Rspec",
                    "Category_Rmax", "pass_worstcase"]].round(3).to_string(index=False))

print("\nCandidates passing the prespecified WORST-CASE rule "
      "(R_max >= 1.20 AND MaxOff < 0.70): %d / %d" % (res.pass_worstcase.sum(), len(res)))
print("Candidates with R_spec > 1.0 (original geometric-mean rule): %d / %d"
      % ((res.R_spec > 1.0).sum(), len(res)))

lead = res[res.scFv == "scfv_73_0"].iloc[0]
print("\n--- Lead candidate scFv-73-0 ---")
print("  on-target CLDN4 STRIVE      : %.3f" % lead.CLDN4)
print("  mean off-target STRIVE      : %.3f  (R_avg = %.3f)" % (lead.AvgOff, lead.R_avg))
print("  worst off-target            : %s = %.3f  (R_max = %.3f)"
      % (lead.MaxOffTarget, lead.MaxOff, lead.R_max))
print("  geometric-mean ratio R_spec : %.3f" % lead.R_spec)
print("  homologues within 10%% of on-target score : %d" % lead.n_off_within_10pct)
print("  homologues failing the 1.20-fold margin   : %d" % lead.n_off_margin_lt_1p20)
print("  passes prespecified worst-case rule       : %s" % lead.pass_worstcase)

print("\n--- scFv-73-0: full ranked off-target profile ---")
prof = mat.loc["scfv_73_0", off_cols].sort_values(ascending=False).to_frame("STRIVE")
prof["margin_vs_CLDN4"] = lead.CLDN4 / prof.STRIVE
print(prof.round(3).to_string())
