#!/usr/bin/env python3
"""
Revision analysis 2 -- Prospective benchmarking of the STRIVE score on an
independent panel of experimentally solved antibody-antigen complexes.
Addresses Reviewer 2 major point 1 (and Reviewer 1 major point 4).

Design (prespecified):
  9 antibody Fv sequences taken from PDB complexes with an experimentally
  determined cognate antigen.  Each Fv is re-predicted with AlphaFold3 (seed 99)
  against its cognate antigen (label = 1) and against 10 non-cognate antigens
  drawn from the other complexes in the panel (label = 0).
  99 antibody-antigen pairs in total; 1 positive and 10 decoys per group.

The STRIVE score is recomputed here exactly as implemented for the manuscript:
  STRIVE = 0.40*ipTM_AbAg + 0.25*phi_PAE + 0.25*phi_affinity + 0.10*pTM - penalty
  phi_PAE      = clip((12 - PAE_interface)/(12 - 5), 0, 1)
  phi_affinity = 0 if dG > -4 else 1/(1 + exp(dG + 9))
  penalty      = max(0, fraction_disordered - 0.15)

Outputs: output/benchmark_per_pair.csv, output/benchmark_per_group.csv
"""
import numpy as np
import pandas as pd

SRC = "/mnt/ssd1/Projects/DesCARTes_wd/scripts/otr_clean_scored_99.csv"
OUT = "output/"
RNG = np.random.default_rng(0)

d = pd.read_csv(SRC)

# --- Recompute STRIVE exactly as in the manuscript implementation -------------
pae_interface = d[["pae_min_HA", "pae_min_LA"]].min(axis=1)   # min over the two Ab chains
phi_pae = np.clip((12.0 - pae_interface) / (12.0 - 5.0), 0, 1)
phi_aff = np.where(d.delta_g > -4.0, 0.0, 1.0 / (1.0 + np.exp(d.delta_g + 9.0)))
penalty = np.maximum(0.0, d.fraction_disordered - 0.15)

d["STRIVE"] = (0.40 * d.iptm_AbAg_mean + 0.25 * phi_pae
               + 0.25 * phi_aff + 0.10 * d.ptmglobal - penalty)
# sensitivity variant using the global ipTM instead of the chain-pair mean
d["STRIVE_globalipTM"] = (0.40 * d.iptm + 0.25 * phi_pae
                          + 0.25 * phi_aff + 0.10 * d.ptmglobal - penalty)
d["phi_PAE"] = phi_pae
d["phi_affinity"] = phi_aff

def auroc(pos, neg):
    """Mann-Whitney U based AUROC."""
    pos, neg = np.asarray(pos), np.asarray(neg)
    allv = np.concatenate([pos, neg])
    r = pd.Series(allv).rank().to_numpy()
    rp = r[: len(pos)].sum()
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))

def boot_auroc_ci(pos, neg, n=10000):
    """Stratified bootstrap CI (resample positives and negatives independently)."""
    pos, neg = np.asarray(pos), np.asarray(neg)
    vals = [auroc(RNG.choice(pos, len(pos), replace=True),
                  RNG.choice(neg, len(neg), replace=True)) for _ in range(n)]
    return np.percentile(vals, [2.5, 97.5])

report = {}
for score in ["STRIVE", "STRIVE_globalipTM", "iptm_AbAg_mean", "phi_PAE",
              "phi_affinity", "ptmglobal"]:
    pos = d.loc[d.label == 1, score].dropna()
    neg = d.loc[d.label == 0, score].dropna()
    a = auroc(pos, neg)
    lo, hi = boot_auroc_ci(pos, neg)
    # per-group rank of the cognate positive (1 = best of 11)
    ranks, top1 = [], 0
    for g, sub in d.groupby("group_id"):
        sub = sub.sort_values(score, ascending=False).reset_index(drop=True)
        rk = int(sub.index[sub.label == 1][0]) + 1
        ranks.append(rk)
        top1 += (rk == 1)
    report[score] = dict(AUROC=a, CI_low=lo, CI_high=hi,
                         top1=top1, n_groups=len(ranks),
                         median_rank=float(np.median(ranks)),
                         mean_reciprocal_rank=float(np.mean([1 / r for r in ranks])),
                         ranks=ranks)

rep = pd.DataFrame(report).T
print("=" * 92)
print("STRIVE BENCHMARK -- 9 cognate antibody-antigen pairs vs 90 non-cognate decoys")
print("=" * 92)
print(rep[["AUROC", "CI_low", "CI_high", "top1", "median_rank",
           "mean_reciprocal_rank"]].astype(float).round(3).to_string())
print("\nRandom expectation: AUROC = 0.500, top-1 = 0.8/9 groups, "
      "median rank = 6, MRR = 0.274")

print("\n--- Per-group rank of the cognate antigen (STRIVE, 1 = best of 11) ---")
per_group = []
for g, sub in d.groupby("group_id"):
    s = sub.sort_values("STRIVE", ascending=False).reset_index(drop=True)
    rk = int(s.index[s.label == 1][0]) + 1
    pos_score = float(s.loc[s.label == 1, "STRIVE"].iloc[0])
    best_decoy = float(s.loc[s.label == 0, "STRIVE"].max())
    per_group.append(dict(group=g, cognate_rank=rk, cognate_STRIVE=pos_score,
                          best_decoy_STRIVE=best_decoy,
                          margin=pos_score - best_decoy,
                          cognate_ranked_first=(rk == 1)))
pg = pd.DataFrame(per_group)
print(pg.round(3).to_string(index=False))

d.to_csv(OUT + "benchmark_per_pair.csv", index=False)
pg.to_csv(OUT + "benchmark_per_group.csv", index=False)
rep.drop(columns=["ranks"]).astype(float).round(4).to_csv(OUT + "benchmark_discrimination.csv")
print("\nSaved -> output/benchmark_{per_pair,per_group,discrimination}.csv")
