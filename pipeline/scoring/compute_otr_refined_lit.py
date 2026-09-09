#!/usr/bin/env python3
"""
compute_otr_refined_lit.py

Add literature-guided antibody off-target metrics to an existing AF3/PRODIGY
OTR results table.

Input:
  - CSV produced by compute_otr_clean.py or af3_otr_pipeline_fixed.py
  - Optional CDR ranges CSV

Main idea:
  Keep AF3 antibody-antigen confidence as the core signal.
  Add penalties for:
    1) antibody sequence/developability liability
    2) PRODIGY/BSA hallucinated geometric fits
    3) poor CDR/paratope-centered engagement

The script does NOT reward high BSA or high raw contact count.
Those are only used as sanity checks / discordance penalties.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


KD_SCALE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

HYDROPHOBIC = set("AVILMFWYC")
AROMATIC = set("FYW")
TYR_TRP = set("YW")
POSITIVE = set("KR")
NEGATIVE = set("DE")


DEFAULT_CDRS = {
    "light": [
        ("L1", 24, 34),
        ("L2", 50, 56),
        ("L3", 89, 97),
    ],
    "heavy": [
        ("H1", 26, 35),
        ("H2", 50, 65),
        ("H3", 95, 102),
    ],
}


def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    if pd.isna(x):
        return np.nan
    return max(lo, min(hi, float(x)))


def sigmoid_affinity(delta_g: float, midpoint: float = -9.0, slope: float = 1.0) -> float:
    if pd.isna(delta_g):
        return np.nan
    return 1.0 / (1.0 + math.exp(slope * (float(delta_g) - midpoint)))


def phi_pae_min(pae_min: float, good: float = 5.0, bad: float = 12.0) -> float:
    """Use pairwise minimum Ab-Ag PAE rather than full mean PAE for this benchmark."""
    if pd.isna(pae_min):
        return np.nan
    pae_min = float(pae_min)
    if pae_min <= good:
        return 1.0
    if pae_min >= bad:
        return 0.0
    return 1.0 - (pae_min - good) / (bad - good)


def infer_group(job_name: str) -> str:
    job = str(job_name).lower()
    if "_lh_vs_" in job:
        return job.split("_lh_vs_")[0]
    if "_vs_" in job:
        return job.split("_vs_")[0]
    return job.split("_")[0]


def infer_label(job_name: str) -> int:
    job = str(job_name).lower()
    if "positive" in job:
        return 1
    if "decoy" in job:
        return 0
    return -1


def find_data_json(job_dir: str | Path) -> Optional[Path]:
    if pd.isna(job_dir):
        return None
    p = Path(str(job_dir))
    if not p.exists():
        return None
    if p.is_file() and p.name.endswith("_data.json"):
        return p
    hits = sorted(p.glob("*_data.json"))
    if hits:
        return hits[0]
    hits = sorted(p.glob("**/*_data.json"))
    if hits:
        return hits[0]
    return None


def load_sequences_from_af3_data(job_dir: str | Path) -> Dict[str, str]:
    """
    Reads AF3 input JSON of the form:
      {"sequences": [{"protein": {"id": "A", "sequence": "..."}}]}
    """
    data_json = find_data_json(job_dir)
    if data_json is None:
        return {}

    try:
        with open(data_json) as f:
            data = json.load(f)
    except Exception:
        return {}

    seqs = {}
    for item in data.get("sequences", []):
        prot = item.get("protein", {})
        cid = prot.get("id")
        seq = prot.get("sequence")
        if isinstance(cid, list):
            # AF3 can support repeated chains. Here we copy same sequence to all IDs.
            for c in cid:
                if seq:
                    seqs[str(c)] = str(seq)
        elif cid is not None and seq:
            seqs[str(cid)] = str(seq)

    return seqs


def read_cdr_ranges(path: Optional[str]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    df = pd.read_csv(path)
    # Accept flexible column names.
    rename = {}
    for c in df.columns:
        lc = c.lower()
        if lc in {"chain_id", "chain"}:
            rename[c] = "chain"
        elif lc in {"cdr", "region", "cdr_name"}:
            rename[c] = "region"
        elif lc in {"start", "start_res", "res_start"}:
            rename[c] = "start"
        elif lc in {"end", "end_res", "res_end"}:
            rename[c] = "end"
    df = df.rename(columns=rename)
    required = {"chain", "region", "start", "end"}
    if not required.issubset(set(df.columns)):
        raise ValueError(
            f"CDR ranges CSV must contain {required}; got {list(df.columns)}"
        )
    if "job_name" not in df.columns:
        df["job_name"] = "*"
    df["job_name"] = df["job_name"].astype(str)
    df["chain"] = df["chain"].astype(str)
    df["region"] = df["region"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def get_ranges_for_job(
    cdr_df: Optional[pd.DataFrame],
    job_name: str,
    light_chain: str,
    heavy_chain: str,
) -> List[Tuple[str, str, int, int]]:
    """
    Returns list of (chain, region, start, end), using 1-based inclusive sequence indices.
    """
    if cdr_df is not None:
        sub = cdr_df[(cdr_df["job_name"] == job_name) | (cdr_df["job_name"] == "*")].copy()
        if len(sub) > 0:
            return [
                (str(r.chain), str(r.region), int(r.start), int(r.end))
                for r in sub.itertuples(index=False)
            ]

    ranges = []
    for region, start, end in DEFAULT_CDRS["light"]:
        ranges.append((light_chain, region, start, end))
    for region, start, end in DEFAULT_CDRS["heavy"]:
        ranges.append((heavy_chain, region, start, end))
    return ranges


def slice_seq(seq: str, start: int, end: int) -> str:
    """1-based inclusive start/end."""
    if not seq:
        return ""
    start0 = max(0, start - 1)
    end0 = min(len(seq), end)
    if start0 >= len(seq) or start0 >= end0:
        return ""
    return seq[start0:end0]


def seq_net_charge(seq: str) -> float:
    return sum(1 for aa in seq if aa in POSITIVE) - sum(1 for aa in seq if aa in NEGATIVE)


def seq_avg_kd(seq: str) -> float:
    vals = [KD_SCALE.get(aa, 0.0) for aa in seq]
    if not vals:
        return np.nan
    return float(np.mean(vals))


def fraction(seq: str, aa_set: set[str]) -> float:
    if not seq:
        return np.nan
    return sum(1 for aa in seq if aa in aa_set) / len(seq)


def compute_sequence_liabilities(row: pd.Series, cdr_df: Optional[pd.DataFrame]) -> Dict[str, float]:
    job_name = str(row.get("job_name", ""))
    job_dir = row.get("job_dir", "")
    light = str(row.get("light_chain", "A"))
    heavy = str(row.get("heavy_chain", "B"))

    seqs = load_sequences_from_af3_data(job_dir)
    ranges = get_ranges_for_job(cdr_df, job_name, light, heavy)

    cdr_records = []
    all_cdr_seq = ""
    hcdr3_seq = ""

    for chain, region, start, end in ranges:
        seq = seqs.get(chain, "")
        frag = slice_seq(seq, start, end)
        all_cdr_seq += frag
        if region.upper() in {"H3", "HCDR3", "CDR-H3"}:
            hcdr3_seq += frag
        cdr_records.append((chain, region, frag))

    cdr_len = len(all_cdr_seq)
    hcdr3_len = len(hcdr3_seq)
    net_charge = seq_net_charge(all_cdr_seq)
    hyd_frac = fraction(all_cdr_seq, HYDROPHOBIC)
    arom_frac = fraction(all_cdr_seq, AROMATIC)
    tyr_trp_frac = fraction(all_cdr_seq, TYR_TRP)
    avg_kd = seq_avg_kd(all_cdr_seq)

    # Literature-guided liability penalties. Conservative by design.
    p_len = 0.00
    if hcdr3_len >= 18:
        p_len = 0.06
    elif hcdr3_len >= 16:
        p_len = 0.04
    elif hcdr3_len >= 14:
        p_len = 0.02

    p_charge = 0.00
    if net_charge > 5:
        p_charge = 0.06
    elif net_charge > 3:
        p_charge = 0.04

    p_hydro = 0.00
    if not pd.isna(hyd_frac):
        if hyd_frac > 0.60:
            p_hydro = 0.06
        elif hyd_frac > 0.45:
            p_hydro = 0.04
        elif hyd_frac > 0.35:
            p_hydro = 0.02

    p_arom = 0.00
    if not pd.isna(tyr_trp_frac):
        if tyr_trp_frac > 0.30:
            p_arom = 0.04
        elif tyr_trp_frac > 0.25:
            p_arom = 0.02

    p_seq = min(0.20, p_len + p_charge + p_hydro + p_arom)

    return {
        "cdr_total_len": cdr_len,
        "hcdr3_len": hcdr3_len,
        "cdr_net_charge": net_charge,
        "cdr_hydrophobic_fraction": hyd_frac,
        "cdr_aromatic_fraction": arom_frac,
        "cdr_tyr_trp_fraction": tyr_trp_frac,
        "cdr_avg_kd": avg_kd,
        "P_seq_liability": p_seq,
        "P_hcdr3_len": p_len,
        "P_cdr_charge": p_charge,
        "P_cdr_hydrophobicity": p_hydro,
        "P_cdr_aromaticity": p_arom,
    }


def compute_discordance_penalties(row: pd.Series) -> Dict[str, float]:
    """
    Penalize cases where geometry/energy looks good but AF3 Ab-Ag confidence is weak.
    """
    iptm_abag = row.get("iptm_AbAg_mean", np.nan)
    pae_min = row.get("pae_min_AbAg_mean", np.nan)
    phi_aff = row.get("phi_affinity_clean", row.get("phi_affinity", np.nan))
    bsa = row.get("bsa_total", np.nan)
    atom_contacts = row.get("atom_contacts_AbAg", np.nan)
    residue_contacts = row.get("residue_contacts_AbAg", np.nan)
    paratope = row.get("paratope_localization", np.nan)

    p_energy_conflict = 0.0
    if not pd.isna(phi_aff):
        if phi_aff > 0.75 and ((not pd.isna(iptm_abag) and iptm_abag < 0.25) or (not pd.isna(pae_min) and pae_min > 12)):
            p_energy_conflict = 0.08
        elif phi_aff > 0.60 and (not pd.isna(iptm_abag) and iptm_abag < 0.20):
            p_energy_conflict = 0.05

    p_bsa_conflict = 0.0
    if not pd.isna(bsa) and not pd.isna(iptm_abag):
        if bsa > 1800 and iptm_abag < 0.25:
            p_bsa_conflict = 0.06
        elif bsa > 2500:
            p_bsa_conflict = 0.04

    p_no_interface = 0.0
    if (not pd.isna(atom_contacts) and atom_contacts < 50) or (not pd.isna(residue_contacts) and residue_contacts < 5):
        p_no_interface = 0.15

    p_paratope = 0.0
    if not pd.isna(paratope):
        # Only meaningful if CDR ranges are proper; harmless if paratope is always 1.0.
        if paratope < 0.50:
            p_paratope = 0.08
        elif paratope < 0.70:
            p_paratope = 0.04

    return {
        "P_energy_conflict": p_energy_conflict,
        "P_bsa_conflict": p_bsa_conflict,
        "P_no_interface_refined": p_no_interface,
        "P_paratope_mislocalization": p_paratope,
    }


def compute_refined_scores(df: pd.DataFrame, cdr_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    df = df.copy()

    if "label" not in df.columns:
        df["label"] = df["job_name"].apply(infer_label)
    if "group_id" not in df.columns:
        df["group_id"] = df["job_name"].apply(infer_group)

    if "phi_affinity_clean" not in df.columns:
        if "delta_g" in df.columns:
            df["phi_affinity_clean"] = df["delta_g"].apply(sigmoid_affinity)
        else:
            df["phi_affinity_clean"] = np.nan

    df["phi_PAE_min_AbAg"] = df["pae_min_AbAg_mean"].apply(phi_pae_min)

    seq_features = []
    discord_features = []
    for _, row in df.iterrows():
        seq_features.append(compute_sequence_liabilities(row, cdr_df))
        discord_features.append(compute_discordance_penalties(row))

    df = pd.concat([df, pd.DataFrame(seq_features), pd.DataFrame(discord_features)], axis=1)

    base_penalty = df["P_penalty_clean"] if "P_penalty_clean" in df.columns else 0.0

    df["P_literature_total"] = (
        df["P_seq_liability"].fillna(0)
        + df["P_energy_conflict"].fillna(0)
        + df["P_bsa_conflict"].fillna(0)
        + df["P_no_interface_refined"].fillna(0)
        + df["P_paratope_mislocalization"].fillna(0)
    )

    # Refined literature-guided OTR:
    # - strong emphasis on antibody-antigen pairwise confidence
    # - PAE uses min Ab-Ag PAE, because mean PAE was saturated at zero in this benchmark
    # - PRODIGY is retained but not allowed to dominate
    # - liability/discordance terms penalize likely nonspecific fits
    df["OTR_refined_lit"] = (
        0.15 * pd.to_numeric(df["iptm"], errors="coerce")
        + 0.40 * pd.to_numeric(df["iptm_AbAg_mean"], errors="coerce")
        + 0.35 * pd.to_numeric(df["phi_PAE_min_AbAg"], errors="coerce")
        + 0.10 * pd.to_numeric(df["phi_affinity_clean"], errors="coerce")
        - base_penalty
        - df["P_literature_total"]
    )

    # No-affinity version for sensitivity analysis.
    df["OTR_refined_lit_noaff"] = (
        0.20 * pd.to_numeric(df["iptm"], errors="coerce")
        + 0.45 * pd.to_numeric(df["iptm_AbAg_mean"], errors="coerce")
        + 0.35 * pd.to_numeric(df["phi_PAE_min_AbAg"], errors="coerce")
        - base_penalty
        - df["P_literature_total"]
    )

    df["rank_OTR_refined_lit"] = (
        df.groupby("group_id")["OTR_refined_lit"]
        .rank(ascending=False, method="min")
        .astype("Int64")
    )
    df["rank_OTR_refined_lit_noaff"] = (
        df.groupby("group_id")["OTR_refined_lit_noaff"]
        .rank(ascending=False, method="min")
        .astype("Int64")
    )

    # Group margins.
    rows = []
    for gid, sub in df.groupby("group_id"):
        pos = sub[sub["label"] == 1]
        dec = sub[sub["label"] == 0]
        pos_score = pos["OTR_refined_lit"].max() if len(pos) else np.nan
        best_dec = dec["OTR_refined_lit"].max() if len(dec) else np.nan
        rows.append({
            "group_id": gid,
            "positive_OTR_refined_lit": pos_score,
            "best_decoy_OTR_refined_lit": best_dec,
            "delta_OTR_refined_lit": pos_score - best_dec if not pd.isna(pos_score) and not pd.isna(best_dec) else np.nan,
        })
    margins = pd.DataFrame(rows)
    df = df.merge(margins, on="group_id", how="left")

    return df


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gid, sub in df.groupby("group_id"):
        pos = sub[sub["label"] == 1]
        if len(pos) == 0:
            continue
        p = pos.sort_values("OTR_refined_lit", ascending=False).iloc[0]
        rows.append({
            "group_id": gid,
            "positive_job": p["job_name"],
            "positive_rank_OTR_refined_lit": p["rank_OTR_refined_lit"],
            "positive_OTR_refined_lit": p["OTR_refined_lit"],
            "best_decoy_OTR_refined_lit": p["best_decoy_OTR_refined_lit"],
            "delta_OTR_refined_lit": p["delta_OTR_refined_lit"],
            "P_seq_liability": p["P_seq_liability"],
            "hcdr3_len": p["hcdr3_len"],
            "cdr_net_charge": p["cdr_net_charge"],
            "cdr_hydrophobic_fraction": p["cdr_hydrophobic_fraction"],
            "cdr_tyr_trp_fraction": p["cdr_tyr_trp_fraction"],
        })
    return pd.DataFrame(rows).sort_values("group_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_csv", required=True, help="Existing scored OTR CSV.")
    ap.add_argument("--output_csv", required=True, help="Output CSV with refined metrics.")
    ap.add_argument("--summary_csv", default=None, help="Optional group-level summary.")
    ap.add_argument("--cdr_ranges", default=None, help="Optional CDR ranges CSV.")
    args = ap.parse_args()

    df = pd.read_csv(args.input_csv)
    cdr_df = read_cdr_ranges(args.cdr_ranges)

    out = compute_refined_scores(df, cdr_df)
    out.to_csv(args.output_csv, index=False)

    print(f"Wrote: {args.output_csv}")
    print(f"Rows: {len(out)}")
    print(f"Rows with OTR_refined_lit: {out['OTR_refined_lit'].notna().sum()}")

    if args.summary_csv:
        sm = summarize(out)
        sm.to_csv(args.summary_csv, index=False)
        print(f"Wrote: {args.summary_csv}")


if __name__ == "__main__":
    main()
