#!/usr/bin/env python3
"""
Compute OTR and antibody-antigen engagement metrics from AlphaFold 3 output directories.

This fixed version supports AF3 output folders like:

job_folder/
  <job_name>_model.cif
  <job_name>_summary_confidences.json
  <job_name>_confidences.json
  ranking_scores.csv
  seed-42_sample-0/
    model.cif
    summary_confidences.json
    confidences.json

Important fixes in this version
-------------------------------
1. --af3_root can be:
   - a single AF3 job directory,
   - a parent directory containing job directories,
   - or a nested tree containing many AF3 job directories.
2. AF3 matrix-style chain_pair_iptm and chain_pair_pae_min are parsed correctly.
   Your AF3 export stores these as NxN arrays, not dictionaries.
3. Mean Ab-Ag PAE is computed from confidences.json using:
   - pae
   - token_chain_ids
4. OTR fallback scores are written even if phi_affinity is unavailable.
5. Optional affinity CSV can be supplied and merged by job_name.

Default chain roles for 3-chain antibody jobs:
  A = light chain
  B = heavy chain
  C = antigen

Typical usage
-------------
Single job parent:
python af3_otr_pipeline_fixed.py \
  --af3_root /path/to/1ahw_LH_vs_1ahw_C_positive/ \
  --output_csv af3_otr_metrics.csv

All jobs recursively:
python af3_otr_pipeline_fixed.py \
  --af3_root /path/to/af3_outputs_99/ \
  --output_csv af3_otr_metrics.csv \
  --recursive

With optional annotations:
python af3_otr_pipeline_fixed.py \
  --af3_root /path/to/af3_outputs_99/ \
  --output_csv af3_otr_metrics.csv \
  --chain_roles chain_roles.csv \
  --cdr_ranges cdr_ranges.csv \
  --affinity_csv affinity.csv \
  --use_freesasa \
  --recursive

Optional chain_roles CSV
------------------------
Columns:
  job_name,light_chain,heavy_chain,antigen_chain

Optional cdr_ranges CSV
-----------------------
Columns:
  job_name,chain_id,start,end,label
Residue positions are 1-based sequence positions within that chain, not PDB numbering.

Optional affinity CSV
---------------------
Columns:
  job_name,phi_affinity
or:
  job_name,affinity_dg
If affinity_dg is supplied, lower/more negative is treated as better and normalized by:
  phi_affinity = clip((bad_dg - affinity_dg) / (bad_dg - good_dg), 0, 1)
Defaults:
  good_dg = -15 kcal/mol
  bad_dg  = 0 kcal/mol
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Any

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

try:
    import freesasa  # type: ignore
    HAVE_FREESASA = True
except Exception:
    HAVE_FREESASA = False


# -----------------------------
# Basic helpers
# -----------------------------

def safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, str) and not x.strip():
            return None
        v = float(x)
        if math.isnan(v):
            return None
        return v
    except Exception:
        return None


def clip01(x: Optional[float]) -> Optional[float]:
    if x is None:
        return None
    return max(0.0, min(1.0, float(x)))


def nanmean(values: Sequence[Optional[float]]) -> Optional[float]:
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if not vals:
        return None
    return float(np.mean(vals))


def norm_higher_better(x: Optional[float], low: float, high: float) -> Optional[float]:
    if x is None or high <= low:
        return None
    return clip01((x - low) / (high - low))


def norm_lower_better(x: Optional[float], good: float, bad: float) -> Optional[float]:
    """1.0 when x <= good; 0.0 when x >= bad."""
    if x is None or bad <= good:
        return None
    return clip01(1.0 - (x - good) / (bad - good))


def choose_first(d: dict, keys: Sequence[str]) -> Optional[float]:
    for k in keys:
        if k in d:
            v = safe_float(d[k])
            if v is not None:
                return v
    return None


def load_json(path: Optional[Path]) -> dict:
    if path is None or not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# -----------------------------
# Configuration dataclasses
# -----------------------------

@dataclass
class ChainRoles:
    light_chain: str = "A"
    heavy_chain: str = "B"
    antigen_chain: str = "C"


@dataclass
class ScoreConfig:
    contact_cutoff: float = 4.5
    proxy_variable_n: int = 120
    pae_good: float = 5.0
    pae_bad: float = 12.0
    contact_low: float = 50.0
    contact_high: float = 800.0
    bsa_low: float = 250.0
    bsa_high: float = 1800.0
    affinity_good_dg: float = -15.0
    affinity_bad_dg: float = 0.0


# -----------------------------
# AF3 job discovery
# -----------------------------

def looks_like_job_dir(path: Path) -> bool:
    if not path.is_dir():
        return False
    has_prefixed = bool(list(path.glob("*_model.cif")) and list(path.glob("*_summary_confidences.json")))
    has_seed_style = (path / "model.cif").exists() and (path / "summary_confidences.json").exists()
    return has_prefixed or has_seed_style


def find_job_dirs(root: Path, recursive: bool = False) -> List[Path]:
    """Return AF3 job directories.

    Accepts:
    - root itself as a job directory
    - direct children as job directories
    - nested trees if recursive=True
    """
    if looks_like_job_dir(root):
        return [root]

    if recursive:
        hits = [p for p in root.rglob("*") if looks_like_job_dir(p)]
    else:
        hits = [p for p in root.iterdir() if p.is_dir() and looks_like_job_dir(p)]

    # Avoid returning both a top-level job and its seed subdirectory if both match.
    hits_sorted = sorted(set(hits), key=lambda p: (len(p.parts), str(p)))
    final: List[Path] = []
    for h in hits_sorted:
        if any(parent in final for parent in h.parents):
            continue
        final.append(h)
    return final


def find_single(job_dir: Path, prefixed_pattern: str, seed_name: str) -> Optional[Path]:
    hits = sorted(job_dir.glob(prefixed_pattern))
    if hits:
        return hits[0]
    seed_path = job_dir / seed_name
    if seed_path.exists():
        return seed_path
    return None


def infer_job_name(job_dir: Path) -> str:
    cifs = sorted(job_dir.glob("*_model.cif"))
    if cifs:
        return re.sub(r"_model$", "", cifs[0].stem)
    data_files = sorted(job_dir.glob("*_data.json"))
    if data_files:
        return re.sub(r"_data$", "", data_files[0].stem)
    return job_dir.name


# -----------------------------
# Chain order and AF3 matrix parsing
# -----------------------------

def unique_in_order(items: Sequence[str]) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()
    for x in items:
        sx = str(x)
        if sx not in seen:
            out.append(sx)
            seen.add(sx)
    return out


def chain_order_from_confidences(conf: dict) -> List[str]:
    token_chains = conf.get("token_chain_ids")
    if isinstance(token_chains, list):
        order = unique_in_order([str(x) for x in token_chains])
        if order:
            return order
    atom_chains = conf.get("atom_chain_ids")
    if isinstance(atom_chains, list):
        order = unique_in_order([str(x) for x in atom_chains])
        if order:
            return order
    return []


def chain_order_from_data(data: dict) -> List[str]:
    # AF3 input data.json usually has sequences: [{protein:{id:'A', sequence:'...'}}]
    seqs = data.get("sequences")
    order: List[str] = []
    if isinstance(seqs, list):
        for rec in seqs:
            if not isinstance(rec, dict):
                continue
            prot = rec.get("protein")
            if isinstance(prot, dict) and "id" in prot:
                cid = prot["id"]
                if isinstance(cid, list):
                    for c in cid:
                        order.append(str(c))
                else:
                    order.append(str(cid))
    return unique_in_order(order)


def infer_chain_order(summary: dict, conf: dict, data: dict, roles: ChainRoles, structure=None) -> List[str]:
    order = chain_order_from_confidences(conf)
    if order:
        return order
    order = chain_order_from_data(data)
    if order:
        return order
    if structure is not None:
        try:
            model = next(structure.get_models())
            order = [chain.id for chain in model]
            if order:
                return order
        except Exception:
            pass
    # Last-resort default matching the generated benchmark JSONs.
    return unique_in_order([roles.light_chain, roles.heavy_chain, roles.antigen_chain])


def matrix_lookup(obj: Any, chain_a: str, chain_b: str, chain_order: Sequence[str]) -> Optional[float]:
    """Look up a pair metric from AF3 summary objects.

    Supports:
    - NxN list matrix indexed by chain_order
    - dict with keys like A_B, A:C
    - nested dict obj[A][C]
    - list of dict records
    """
    if obj is None:
        return None

    # AF3 current export: list-of-lists matrix.
    if isinstance(obj, list) and obj and all(isinstance(r, list) for r in obj):
        idx = {cid: i for i, cid in enumerate(chain_order)}
        if chain_a in idx and chain_b in idx:
            i, j = idx[chain_a], idx[chain_b]
            try:
                return safe_float(obj[i][j])
            except Exception:
                return None

    # Flat or nested dictionaries.
    if isinstance(obj, dict):
        for k in [
            f"{chain_a}_{chain_b}", f"{chain_b}_{chain_a}",
            f"{chain_a}:{chain_b}", f"{chain_b}:{chain_a}",
            f"{chain_a}-{chain_b}", f"{chain_b}-{chain_a}",
        ]:
            if k in obj:
                return safe_float(obj[k])
        if chain_a in obj and isinstance(obj[chain_a], dict) and chain_b in obj[chain_a]:
            return safe_float(obj[chain_a][chain_b])
        if chain_b in obj and isinstance(obj[chain_b], dict) and chain_a in obj[chain_b]:
            return safe_float(obj[chain_b][chain_a])

    # List of records.
    if isinstance(obj, list):
        for rec in obj:
            if not isinstance(rec, dict):
                continue
            vals = {str(v) for v in rec.values() if isinstance(v, (str, int))}
            if chain_a in vals and chain_b in vals:
                for metric_key in ["value", "iptm", "ipTM", "pae_min", "paeMin", "pae_mean", "paeMean", "score"]:
                    if metric_key in rec:
                        v = safe_float(rec[metric_key])
                        if v is not None:
                            return v
    return None


def extract_global_summary_metrics(summary: dict) -> Dict[str, Optional[float]]:
    return {
        "iptm": choose_first(summary, ["iptm", "ipTM"]),
        "ptmglobal": choose_first(summary, ["ptm", "pTM"]),
        "ranking_score": choose_first(summary, ["ranking_score", "rankingScore"]),
        "fraction_disordered": choose_first(summary, ["fraction_disordered", "fractionDisordered"]),
        "has_clash": choose_first(summary, ["has_clash", "hasClash"]),
    }


def extract_pairwise_summary_metrics(summary: dict, roles: ChainRoles, chain_order: Sequence[str]) -> Dict[str, Optional[float]]:
    pair_iptm_obj = summary.get("chain_pair_iptm") or summary.get("chainPairIptm")
    pair_pae_min_obj = summary.get("chain_pair_pae_min") or summary.get("chainPairPaeMin")
    pair_pae_mean_obj = summary.get("chain_pair_pae_mean") or summary.get("chainPairPaeMean")

    out: Dict[str, Optional[float]] = {
        "iptm_HA": matrix_lookup(pair_iptm_obj, roles.heavy_chain, roles.antigen_chain, chain_order),
        "iptm_LA": matrix_lookup(pair_iptm_obj, roles.light_chain, roles.antigen_chain, chain_order),
        "pae_min_HA": matrix_lookup(pair_pae_min_obj, roles.heavy_chain, roles.antigen_chain, chain_order),
        "pae_min_LA": matrix_lookup(pair_pae_min_obj, roles.light_chain, roles.antigen_chain, chain_order),
        "pae_mean_HA": matrix_lookup(pair_pae_mean_obj, roles.heavy_chain, roles.antigen_chain, chain_order),
        "pae_mean_LA": matrix_lookup(pair_pae_mean_obj, roles.light_chain, roles.antigen_chain, chain_order),
    }
    out["iptm_AbAg_mean"] = nanmean([out["iptm_HA"], out["iptm_LA"]])
    out["pae_min_AbAg_mean"] = nanmean([out["pae_min_HA"], out["pae_min_LA"]])
    out["pae_mean_AbAg_mean"] = nanmean([out["pae_mean_HA"], out["pae_mean_LA"]])
    return out


def mean_pae_between_chains(conf: dict, chain1: str, chain2: str) -> Optional[float]:
    pae = conf.get("pae")
    token_chain_ids = conf.get("token_chain_ids")
    if pae is None or token_chain_ids is None:
        return None
    if not isinstance(token_chain_ids, list):
        return None

    idx1 = [i for i, c in enumerate(token_chain_ids) if str(c) == chain1]
    idx2 = [i for i, c in enumerate(token_chain_ids) if str(c) == chain2]
    if not idx1 or not idx2:
        return None

    # Avoid materializing huge submatrices for very large jobs; loop is fine for typical Fab/antigen.
    vals: List[float] = []
    try:
        for i in idx1:
            row_i = pae[i]
            for j in idx2:
                v1 = safe_float(row_i[j])
                if v1 is not None:
                    vals.append(v1)
                v2 = safe_float(pae[j][i])
                if v2 is not None:
                    vals.append(v2)
    except Exception:
        return None

    if not vals:
        return None
    return float(np.mean(vals))


def extract_mean_pae_metrics(conf: dict, roles: ChainRoles) -> Dict[str, Optional[float]]:
    pae_mean_LA = mean_pae_between_chains(conf, roles.light_chain, roles.antigen_chain)
    pae_mean_HA = mean_pae_between_chains(conf, roles.heavy_chain, roles.antigen_chain)
    return {
        "pae_mean_LA": pae_mean_LA,
        "pae_mean_HA": pae_mean_HA,
        "pae_mean_AbAg_mean": nanmean([pae_mean_HA, pae_mean_LA]),
    }


# -----------------------------
# Structure parsing and contacts
# -----------------------------

def parse_structure(cif_path: Path):
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(cif_path.stem, str(cif_path))


def get_chain(structure, chain_id: str):
    model = next(structure.get_models())
    if chain_id not in model:
        available = [c.id for c in model]
        raise KeyError(f"Chain {chain_id!r} not found. Available chains: {available}")
    return model[chain_id]


def is_protein_residue(res: Residue) -> bool:
    return res.id[0] == " "


def residue_key(res: Residue) -> Tuple[str, int, str]:
    return (res.get_parent().id, int(res.id[1]), str(res.id[2]).strip() or "")


def iter_heavy_atoms(chain, include_hydrogens: bool = False) -> Iterable[Atom]:
    for res in chain:
        if not is_protein_residue(res):
            continue
        for atom in res.get_atoms():
            element = (atom.element or atom.get_name()[0]).strip().upper()
            if not include_hydrogens and element == "H":
                continue
            yield atom


def unique_residue_contacts(chain_a, chain_b, cutoff: float) -> Tuple[Set[Tuple], Set[Tuple], int]:
    atoms_b = list(iter_heavy_atoms(chain_b))
    ns = NeighborSearch(atoms_b)
    residues_a: Set[Tuple] = set()
    residues_b: Set[Tuple] = set()
    atom_pairs = 0

    for atom_a in iter_heavy_atoms(chain_a):
        nearby = ns.search(atom_a.coord, cutoff, level="A")
        for atom_b in nearby:
            atom_pairs += 1
            residues_a.add(residue_key(atom_a.get_parent()))
            residues_b.add(residue_key(atom_b.get_parent()))

    return residues_a, residues_b, atom_pairs


def chain_length(chain) -> int:
    return sum(1 for res in chain if is_protein_residue(res))


def residue_index_map(chain) -> List[Tuple[Tuple[str, int, str], int]]:
    out: List[Tuple[Tuple[str, int, str], int]] = []
    idx = 0
    for res in chain:
        if not is_protein_residue(res):
            continue
        idx += 1
        out.append((residue_key(res), idx))
    return out


def proxy_variable_domain_residues(chain, n_res: int) -> Set[Tuple]:
    return {rk for rk, idx in residue_index_map(chain) if idx <= n_res}


# -----------------------------
# Optional CDR ranges and roles
# -----------------------------

def load_chain_roles_csv(path: Optional[Path]) -> Dict[str, ChainRoles]:
    if path is None:
        return {}
    df = pd.read_csv(path)
    required = {"job_name", "light_chain", "heavy_chain", "antigen_chain"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"chain_roles CSV missing columns: {sorted(missing)}")
    out: Dict[str, ChainRoles] = {}
    for _, row in df.iterrows():
        out[str(row["job_name"])] = ChainRoles(
            light_chain=str(row["light_chain"]),
            heavy_chain=str(row["heavy_chain"]),
            antigen_chain=str(row["antigen_chain"]),
        )
    return out


def load_cdr_ranges_csv(path: Optional[Path]) -> Dict[str, Dict[str, List[Tuple[int, int, str]]]]:
    if path is None:
        return {}
    df = pd.read_csv(path)
    required = {"job_name", "chain_id", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"cdr_ranges CSV missing columns: {sorted(missing)}")
    out: Dict[str, Dict[str, List[Tuple[int, int, str]]]] = defaultdict(lambda: defaultdict(list))
    for _, row in df.iterrows():
        label = str(row["label"]) if "label" in df.columns else ""
        out[str(row["job_name"])][str(row["chain_id"])].append((int(row["start"]), int(row["end"]), label))
    return out


def cdr_residue_set(chain, ranges: List[Tuple[int, int, str]]) -> Set[Tuple]:
    seq_to_key = {idx: rk for rk, idx in residue_index_map(chain)}
    out: Set[Tuple] = set()
    for start, end, _label in ranges:
        for idx in range(start, end + 1):
            rk = seq_to_key.get(idx)
            if rk is not None:
                out.add(rk)
    return out


# -----------------------------
# Optional affinity input
# -----------------------------

def load_affinity_csv(path: Optional[Path], cfg: ScoreConfig) -> Dict[str, float]:
    if path is None:
        return {}
    df = pd.read_csv(path)
    if "job_name" not in df.columns:
        raise ValueError("affinity CSV must contain job_name")
    if "phi_affinity" not in df.columns and "affinity_dg" not in df.columns:
        raise ValueError("affinity CSV must contain phi_affinity or affinity_dg")
    out: Dict[str, float] = {}
    for _, row in df.iterrows():
        job = str(row["job_name"])
        if "phi_affinity" in df.columns and pd.notna(row.get("phi_affinity")):
            val = clip01(safe_float(row["phi_affinity"]))
        else:
            dg = safe_float(row.get("affinity_dg"))
            val = None
            if dg is not None:
                # More negative is better. good_dg gets 1; bad_dg gets 0.
                val = clip01((cfg.affinity_bad_dg - dg) / (cfg.affinity_bad_dg - cfg.affinity_good_dg))
        if val is not None:
            out[job] = val
    return out


# -----------------------------
# Optional FreeSASA BSA
# -----------------------------

def write_chain_subset_pdb(structure, chain_ids: Sequence[str], out_path: Path) -> None:
    model = next(structure.get_models())
    serial = 1
    with open(out_path, "w", encoding="utf-8") as f:
        for chain_id in chain_ids:
            if chain_id not in model:
                continue
            chain = model[chain_id]
            for res in chain:
                if not is_protein_residue(res):
                    continue
                for atom in res.get_atoms():
                    element = (atom.element or atom.get_name()[0]).strip().upper()
                    if element == "H":
                        continue
                    name = atom.get_name()
                    altloc = atom.get_altloc() if atom.get_altloc() != " " else ""
                    x, y, z = atom.coord
                    resname = res.resname
                    chain_char = chain.id[:1]
                    resseq = int(res.id[1])
                    icode = (res.id[2] or " ").strip() or " "
                    occ = atom.occupancy if atom.occupancy is not None else 1.0
                    bfac = atom.bfactor if atom.bfactor is not None else 0.0
                    f.write(
                        f"ATOM  {serial:5d} {name:<4}{altloc:1}{resname:>3} {chain_char:1}"
                        f"{resseq:4d}{icode:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          {element:>2}\n"
                    )
                    serial += 1
        f.write("END\n")


def compute_sasa(path: Path) -> float:
    st = freesasa.Structure(str(path))
    result = freesasa.calc(st)
    return float(result.totalArea())


def compute_bsa(structure, roles: ChainRoles, tmp_dir: Path) -> Dict[str, Optional[float]]:
    if not HAVE_FREESASA:
        return {"bsa_total": None, "bsa_HA": None, "bsa_LA": None}

    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "complex": tmp_dir / "complex.pdb",
        "ab": tmp_dir / "ab.pdb",
        "ag": tmp_dir / "ag.pdb",
        "ha": tmp_dir / "ha.pdb",
        "la": tmp_dir / "la.pdb",
        "h": tmp_dir / "h.pdb",
        "l": tmp_dir / "l.pdb",
    }
    write_chain_subset_pdb(structure, [roles.light_chain, roles.heavy_chain, roles.antigen_chain], paths["complex"])
    write_chain_subset_pdb(structure, [roles.light_chain, roles.heavy_chain], paths["ab"])
    write_chain_subset_pdb(structure, [roles.antigen_chain], paths["ag"])
    write_chain_subset_pdb(structure, [roles.heavy_chain, roles.antigen_chain], paths["ha"])
    write_chain_subset_pdb(structure, [roles.light_chain, roles.antigen_chain], paths["la"])
    write_chain_subset_pdb(structure, [roles.heavy_chain], paths["h"])
    write_chain_subset_pdb(structure, [roles.light_chain], paths["l"])

    try:
        sasa_ab = compute_sasa(paths["ab"])
        sasa_ag = compute_sasa(paths["ag"])
        sasa_complex = compute_sasa(paths["complex"])
        sasa_h = compute_sasa(paths["h"])
        sasa_l = compute_sasa(paths["l"])
        sasa_ha = compute_sasa(paths["ha"])
        sasa_la = compute_sasa(paths["la"])
        return {
            "bsa_total": sasa_ab + sasa_ag - sasa_complex,
            "bsa_HA": sasa_h + sasa_ag - sasa_ha,
            "bsa_LA": sasa_l + sasa_ag - sasa_la,
        }
    except Exception:
        return {"bsa_total": None, "bsa_HA": None, "bsa_LA": None}


# -----------------------------
# Main per-job metric calculation
# -----------------------------

def compute_job_metrics(
    job_dir: Path,
    roles: ChainRoles,
    cdr_config: Dict[str, Dict[str, List[Tuple[int, int, str]]]],
    affinity_map: Dict[str, float],
    use_freesasa: bool,
    cfg: ScoreConfig,
) -> Dict[str, Any]:
    job_name = infer_job_name(job_dir)

    cif_path = find_single(job_dir, "*_model.cif", "model.cif")
    summary_path = find_single(job_dir, "*_summary_confidences.json", "summary_confidences.json")
    conf_path = find_single(job_dir, "*_confidences.json", "confidences.json")
    data_path = find_single(job_dir, "*_data.json", "data.json")

    if cif_path is None or summary_path is None:
        raise FileNotFoundError(f"Missing model.cif or summary_confidences.json in {job_dir}")

    summary = load_json(summary_path)
    conf = load_json(conf_path)
    data = load_json(data_path)

    structure = parse_structure(cif_path)
    chain_order = infer_chain_order(summary, conf, data, roles, structure=structure)

    light = get_chain(structure, roles.light_chain)
    heavy = get_chain(structure, roles.heavy_chain)
    antigen = get_chain(structure, roles.antigen_chain)

    # AF3 confidence metrics.
    base = extract_global_summary_metrics(summary)
    pair = extract_pairwise_summary_metrics(summary, roles, chain_order)

    # Full PAE mean from confidences.json overrides missing summary mean values.
    mean_pae = extract_mean_pae_metrics(conf, roles)
    for k, v in mean_pae.items():
        if v is not None:
            pair[k] = v

    # Geometry contacts.
    l_contact_res, ag_res_from_l, atom_contacts_LA = unique_residue_contacts(light, antigen, cfg.contact_cutoff)
    h_contact_res, ag_res_from_h, atom_contacts_HA = unique_residue_contacts(heavy, antigen, cfg.contact_cutoff)
    ag_contact_res = ag_res_from_l | ag_res_from_h

    residue_contacts_LA = len(l_contact_res)
    residue_contacts_HA = len(h_contact_res)
    residue_contacts_AbAg = residue_contacts_LA + residue_contacts_HA
    atom_contacts_AbAg = atom_contacts_LA + atom_contacts_HA

    antigen_len = chain_length(antigen)
    antigen_contact_fraction = len(ag_contact_res) / antigen_len if antigen_len > 0 else None

    denom = atom_contacts_HA + atom_contacts_LA
    contact_balance_HL = 1.0 - abs(atom_contacts_HA - atom_contacts_LA) / denom if denom > 0 else None

    # Paratope localization.
    if job_name in cdr_config:
        l_para = cdr_residue_set(light, cdr_config[job_name].get(roles.light_chain, []))
        h_para = cdr_residue_set(heavy, cdr_config[job_name].get(roles.heavy_chain, []))
        paratope_mode = "explicit_cdr"
    else:
        l_para = proxy_variable_domain_residues(light, cfg.proxy_variable_n)
        h_para = proxy_variable_domain_residues(heavy, cfg.proxy_variable_n)
        paratope_mode = f"proxy_first_{cfg.proxy_variable_n}res"

    localized_contacts = len(l_contact_res & l_para) + len(h_contact_res & h_para)
    paratope_localization = localized_contacts / residue_contacts_AbAg if residue_contacts_AbAg > 0 else None

    # Optional BSA.
    if use_freesasa:
        bsa = compute_bsa(structure, roles, job_dir / ".tmp_freesasa")
    else:
        bsa = {"bsa_total": None, "bsa_HA": None, "bsa_LA": None}

    # Normalized components.
    pae_for_phi = pair.get("pae_mean_AbAg_mean")
    if pae_for_phi is None:
        pae_for_phi = pair.get("pae_min_AbAg_mean")

    phi_PAE = norm_lower_better(pae_for_phi, good=cfg.pae_good, bad=cfg.pae_bad)
    phi_AbAg_PAE = phi_PAE
    phi_contacts = norm_higher_better(atom_contacts_AbAg, cfg.contact_low, cfg.contact_high)
    phi_bsa = norm_higher_better(bsa.get("bsa_total"), cfg.bsa_low, cfg.bsa_high)
    phi_balance = clip01(contact_balance_HL)
    phi_paratope = clip01(paratope_localization)

    # Weighted engagement correction. Missing terms are not imputed; denominator is renormalized.
    engagement_terms = [
        (0.35, phi_AbAg_PAE),
        (0.25, phi_contacts),
        (0.20, phi_bsa),
        (0.10, phi_balance),
        (0.10, phi_paratope),
    ]
    num = sum(w * v for w, v in engagement_terms if v is not None)
    den = sum(w for w, v in engagement_terms if v is not None)
    phi_engagement = num / den if den > 0 else None

    # Affinity can be merged from external CSV. If absent, keep phi_affinity blank but compute no-affinity scores.
    phi_affinity = affinity_map.get(job_name)

    # Penalty hooks.
    penalty = 0.0
    if base.get("has_clash") is not None and float(base["has_clash"]) > 0:
        penalty += 0.05
    if base.get("fraction_disordered") is not None and float(base["fraction_disordered"]) > 0.30:
        penalty += 0.05

    iptm = base.get("iptm")
    ptm = base.get("ptmglobal")

    # Original OTR with affinity if available, and no-affinity fallback.
    otr_original = None
    otr_original_noaff = None
    if iptm is not None and phi_PAE is not None and ptm is not None:
        otr_original_noaff = 0.40 * iptm + 0.25 * phi_PAE + 0.10 * ptm - penalty
        aff = 0.0 if phi_affinity is None else phi_affinity
        otr_original = 0.40 * iptm + 0.25 * phi_PAE + 0.25 * aff + 0.10 * ptm - penalty

    # Fab-aware version; these are proposed weights for calibration, not final fitted weights.
    otr_fab = None
    otr_fab_noaff = None
    if iptm is not None and phi_PAE is not None and ptm is not None and phi_engagement is not None:
        otr_fab_noaff = 0.30 * iptm + 0.20 * phi_PAE + 0.05 * ptm + 0.35 * phi_engagement - penalty
        aff = 0.0 if phi_affinity is None else phi_affinity
        otr_fab = 0.30 * iptm + 0.20 * phi_PAE + 0.10 * aff + 0.05 * ptm + 0.35 * phi_engagement - penalty

    return {
        "job_name": job_name,
        "job_dir": str(job_dir),
        "model_cif": str(cif_path),
        "summary_json": str(summary_path),
        "confidences_json": str(conf_path) if conf_path else "",
        "chain_order": ";".join(chain_order),
        "light_chain": roles.light_chain,
        "heavy_chain": roles.heavy_chain,
        "antigen_chain": roles.antigen_chain,
        "paratope_mode": paratope_mode,
        **base,
        **pair,
        "atom_contacts_LA": atom_contacts_LA,
        "atom_contacts_HA": atom_contacts_HA,
        "atom_contacts_AbAg": atom_contacts_AbAg,
        "residue_contacts_LA": residue_contacts_LA,
        "residue_contacts_HA": residue_contacts_HA,
        "residue_contacts_AbAg": residue_contacts_AbAg,
        "antigen_residues_contacted": len(ag_contact_res),
        "antigen_length": antigen_len,
        "antigen_contact_fraction": antigen_contact_fraction,
        "contact_balance_HL": contact_balance_HL,
        "paratope_localization": paratope_localization,
        "bsa_total": bsa.get("bsa_total"),
        "bsa_HA": bsa.get("bsa_HA"),
        "bsa_LA": bsa.get("bsa_LA"),
        "phi_PAE": phi_PAE,
        "phi_affinity": phi_affinity,
        "penalty": penalty,
        "phi_AbAg_PAE": phi_AbAg_PAE,
        "phi_contacts": phi_contacts,
        "phi_bsa": phi_bsa,
        "phi_balance": phi_balance,
        "phi_paratope": phi_paratope,
        "phi_engagement": phi_engagement,
        "otr_original": otr_original,
        "otr_original_noaff": otr_original_noaff,
        "otr_fab": otr_fab,
        "otr_fab_noaff": otr_fab_noaff,
    }


# -----------------------------
# CLI
# -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute OTR and Fab engagement metrics from AF3 outputs.")
    p.add_argument("--af3_root", required=True, help="AF3 job directory, parent directory, or root tree.")
    p.add_argument("--output_csv", required=True, help="Output CSV path.")
    p.add_argument("--recursive", action="store_true", help="Search recursively for AF3 job directories under --af3_root.")
    p.add_argument("--chain_roles", default=None, help="Optional CSV mapping job_name to light/heavy/antigen chains.")
    p.add_argument("--cdr_ranges", default=None, help="Optional CSV of CDR/paratope residue ranges.")
    p.add_argument("--affinity_csv", default=None, help="Optional CSV with job_name plus phi_affinity or affinity_dg.")
    p.add_argument("--contact_cutoff", type=float, default=4.5, help="Heavy-atom contact cutoff in Angstrom.")
    p.add_argument("--proxy_variable_n", type=int, default=120, help="If no CDR CSV, first N H/L residues are paratope proxy.")
    p.add_argument("--pae_good", type=float, default=5.0, help="PAE value receiving phi_PAE=1.")
    p.add_argument("--pae_bad", type=float, default=12.0, help="PAE value receiving phi_PAE=0.")
    p.add_argument("--use_freesasa", action="store_true", help="Compute BSA using FreeSASA if installed.")
    p.add_argument("--affinity_good_dg", type=float, default=-15.0, help="Affinity dG mapped to phi_affinity=1.")
    p.add_argument("--affinity_bad_dg", type=float, default=0.0, help="Affinity dG mapped to phi_affinity=0.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.af3_root).expanduser().resolve()
    output_csv = Path(args.output_csv).expanduser().resolve()

    if not root.exists() or not root.is_dir():
        print(f"ERROR: --af3_root is not a directory: {root}", file=sys.stderr)
        return 2

    cfg = ScoreConfig(
        contact_cutoff=args.contact_cutoff,
        proxy_variable_n=args.proxy_variable_n,
        pae_good=args.pae_good,
        pae_bad=args.pae_bad,
        affinity_good_dg=args.affinity_good_dg,
        affinity_bad_dg=args.affinity_bad_dg,
    )

    if args.use_freesasa and not HAVE_FREESASA:
        print("[WARN] --use_freesasa requested but freesasa is not installed. BSA columns will be blank.", file=sys.stderr)

    roles_map = load_chain_roles_csv(Path(args.chain_roles).expanduser().resolve() if args.chain_roles else None)
    cdr_cfg = load_cdr_ranges_csv(Path(args.cdr_ranges).expanduser().resolve() if args.cdr_ranges else None)
    affinity_map = load_affinity_csv(Path(args.affinity_csv).expanduser().resolve() if args.affinity_csv else None, cfg)

    job_dirs = find_job_dirs(root, recursive=args.recursive)
    if not job_dirs:
        print(f"ERROR: no AF3 job directories found under {root}", file=sys.stderr)
        return 2

    rows: List[Dict[str, Any]] = []
    for job_dir in job_dirs:
        job_name = infer_job_name(job_dir)
        roles = roles_map.get(job_name, ChainRoles())
        try:
            row = compute_job_metrics(
                job_dir=job_dir,
                roles=roles,
                cdr_config=cdr_cfg,
                affinity_map=affinity_map,
                use_freesasa=args.use_freesasa,
                cfg=cfg,
            )
            rows.append(row)
            print(f"[OK] {job_name}")
        except Exception as exc:
            print(f"[WARN] {job_name}: {exc}", file=sys.stderr)

    if not rows:
        print("ERROR: no jobs processed successfully", file=sys.stderr)
        return 2

    df = pd.DataFrame(rows)
    preferred = [
        "job_name", "job_dir", "model_cif", "summary_json", "confidences_json",
        "chain_order", "light_chain", "heavy_chain", "antigen_chain", "paratope_mode",
        "iptm", "ptmglobal", "ranking_score", "fraction_disordered", "has_clash",
        "iptm_HA", "iptm_LA", "iptm_AbAg_mean",
        "pae_min_HA", "pae_min_LA", "pae_min_AbAg_mean",
        "pae_mean_HA", "pae_mean_LA", "pae_mean_AbAg_mean",
        "atom_contacts_LA", "atom_contacts_HA", "atom_contacts_AbAg",
        "residue_contacts_LA", "residue_contacts_HA", "residue_contacts_AbAg",
        "antigen_residues_contacted", "antigen_length", "antigen_contact_fraction",
        "contact_balance_HL", "paratope_localization",
        "bsa_total", "bsa_HA", "bsa_LA",
        "phi_PAE", "phi_affinity", "penalty",
        "phi_AbAg_PAE", "phi_contacts", "phi_bsa", "phi_balance", "phi_paratope", "phi_engagement",
        "otr_original", "otr_original_noaff", "otr_fab", "otr_fab_noaff",
    ]
    ordered_cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    df = df[ordered_cols]

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"[DONE] wrote {output_csv} with {len(df)} rows")

    readme_path = output_csv.with_suffix(".README.txt")
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(
            "AF3 OTR pipeline output notes\n"
            "=============================\n\n"
            "Key fixes in this version:\n"
            "- chain_pair_iptm and chain_pair_pae_min matrix outputs are indexed by chain_order.\n"
            "- chain_order is inferred from token_chain_ids in confidences.json where available.\n"
            "- pae_mean_* is computed from the full token-level PAE matrix.\n"
            "- otr_original_noaff and otr_fab_noaff are provided when no phi_affinity is supplied.\n\n"
            "Recommended score columns while affinity is absent:\n"
            "- otr_original_noaff for original OTR without affinity.\n"
            "- otr_fab_noaff for Fab-aware OTR without affinity.\n\n"
            "If you later add affinity values, provide --affinity_csv with job_name,phi_affinity.\n"
        )
    print(f"[DONE] wrote {readme_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
