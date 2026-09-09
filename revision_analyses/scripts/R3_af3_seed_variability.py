#!/usr/bin/env python3
"""
Revision analysis 3 -- AlphaFold3 run-to-run (seed-to-seed) variability of the
STRIVE structural terms for antibody-claudin complexes.
Addresses Reviewer 2 major points 2 and 6.

Five scFv-claudin pairs were predicted twice with independent AlphaFold3 seeds
(seed 10 and seed 42).  We report the seed-to-seed difference in the structural
components of STRIVE, which bounds the resolution at which two STRIVE scores can
be meaningfully distinguished.

Outputs: output/af3_seed_variability.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/mnt/ssd1/Projects/DesCARTes_wd/binder/RFantibody/CLDN_family/af3_test_output")
OUT = "output/"

def parse(p):
    d = json.load(open(p))
    cp_iptm = d.get("chain_pair_iptm", [])
    cp_pae = d.get("chain_pair_pae_min", [])
    iptm = (cp_iptm[0][2] + cp_iptm[1][2]) / 2 if len(cp_iptm) >= 3 else d.get("iptm", np.nan)
    pae = min(cp_pae[0][2], cp_pae[1][2]) if len(cp_pae) >= 3 else np.nan
    return dict(iptm_AbAg=iptm, pae_interface=pae, ptm=d.get("ptm", np.nan),
                fraction_disordered=d.get("fraction_disordered", 0.0),
                has_clash=d.get("has_clash", 0))

rows = []
for job in sorted(ROOT.iterdir()):
    if not job.is_dir():
        continue
    for seed_dir in sorted(job.glob("seed-*_sample-*")):
        f = seed_dir / "summary_confidences.json"
        if not f.exists():
            continue
        seed = int(seed_dir.name.split("_")[0].split("-")[1])
        rows.append(dict(pair=job.name, seed=seed, **parse(f)))

d = pd.DataFrame(rows)
# structural part of STRIVE (affinity term excluded: PRODIGY dG not recomputed here)
d["phi_PAE"] = np.clip((12.0 - d.pae_interface) / 7.0, 0, 1)
d["penalty"] = np.maximum(0.0, d.fraction_disordered - 0.15)
d["STRIVE_struct"] = (0.40 * d.iptm_AbAg + 0.25 * d.phi_PAE
                      + 0.10 * d.ptm - d.penalty)
# renormalised to the 0.75 of total weight that is structural, so the numbers are
# on the same 0-1 scale as the published STRIVE score
d["STRIVE_struct_scaled"] = d.STRIVE_struct / 0.75

w = d.pivot(index="pair", columns="seed",
            values=["iptm_AbAg", "pae_interface", "ptm", "STRIVE_struct_scaled"])
w.columns = [f"{a}_seed{b}" for a, b in w.columns]
for m in ["iptm_AbAg", "pae_interface", "ptm", "STRIVE_struct_scaled"]:
    w[f"{m}_delta"] = (w[f"{m}_seed42"] - w[f"{m}_seed10"]).abs()

pd.set_option("display.width", 200)
print("=" * 88)
print("AlphaFold3 SEED-TO-SEED VARIABILITY  (n = 5 scFv-claudin pairs, seeds 10 vs 42)")
print("=" * 88)
print(w.round(3).to_string())

print("\n--- Absolute seed-to-seed differences ---")
for m, lab in [("iptm_AbAg", "interface ipTM"),
               ("pae_interface", "interface PAE (A)"),
               ("ptm", "global pTM"),
               ("STRIVE_struct_scaled", "STRIVE (structural terms, rescaled)")]:
    v = w[f"{m}_delta"]
    print(f"  {lab:38s} mean |delta| = {v.mean():.3f}   max |delta| = {v.max():.3f}")

sd = w["STRIVE_struct_scaled_delta"]
print("\nInterpretation: with a mean absolute seed-to-seed shift of "
      f"{sd.mean():.3f} STRIVE units (max {sd.max():.3f}), a difference between two "
      "STRIVE scores")
print("smaller than ~%.2f units cannot be distinguished from AlphaFold3 sampling noise."
      % sd.max())
print("The scFv-73-0 CLDN4 (0.846) vs CLDN18.2 (0.798) difference is %.3f units."
      % (0.846 - 0.798))

w.round(4).to_csv(OUT + "af3_seed_variability.csv")
d.round(4).to_csv(OUT + "af3_seed_variability_long.csv", index=False)
print("\nSaved -> output/af3_seed_variability.csv")
