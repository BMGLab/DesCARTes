#!/usr/bin/env python3

import argparse
import math
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


def sigmoid_affinity(delta_g: float, midpoint: float = -9.0, slope: float = 1.0) -> float:
    if pd.isna(delta_g):
        return float("nan")
    return 1.0 / (1.0 + math.exp(slope * (delta_g - midpoint)))


def phi_pae(pae_value: float, good: float = 5.0, bad: float = 12.0) -> float:
    if pd.isna(pae_value):
        return float("nan")
    if pae_value <= good:
        return 1.0
    if pae_value >= bad:
        return 0.0
    return 1.0 - ((pae_value - good) / (bad - good))


def choose_pae(row: pd.Series) -> float:
    for col in ["pae_mean_AbAg_mean", "pae_min_AbAg_mean"]:
        if col in row and not pd.isna(row[col]):
            return row[col]
    return float("nan")


def compute_penalty(
    row: pd.Series,
    disorder_threshold: float = 0.15,
    min_atom_contacts: int = 50,
    min_residue_contacts: int = 5,
    use_bsa_penalty: bool = False,
    min_bsa: float = 400.0,
) -> float:
    penalty = 0.0

    fraction_disordered = row.get("fraction_disordered", float("nan"))
    if not pd.isna(fraction_disordered) and fraction_disordered > disorder_threshold:
        penalty += 0.10

    has_clash = row.get("has_clash", 0)
    if not pd.isna(has_clash) and float(has_clash) > 0:
        penalty += 0.15

    atom_contacts = row.get("atom_contacts_AbAg", float("nan"))
    residue_contacts = row.get("residue_contacts_AbAg", float("nan"))

    no_interface = False
    if not pd.isna(atom_contacts) and atom_contacts < min_atom_contacts:
        no_interface = True
    if not pd.isna(residue_contacts) and residue_contacts < min_residue_contacts:
        no_interface = True

    if no_interface:
        penalty += 0.15

    if use_bsa_penalty and "bsa_total" in row:
        bsa = row.get("bsa_total", float("nan"))
        if not pd.isna(bsa) and bsa < min_bsa:
            penalty += 0.10

    return penalty


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


def load_affinity_table(path: str | None) -> pd.DataFrame | None:
    if path is None:
        return None

    aff = pd.read_csv(path)

    if "job_name" not in aff.columns:
        raise ValueError("Affinity CSV must contain a 'job_name' column.")

    dg_col = None
    for candidate in ["delta_g", "prodigy_delta_g", "DeltaG", "dG", "dg"]:
        if candidate in aff.columns:
            dg_col = candidate
            break

    if dg_col is None:
        raise ValueError(
            "Affinity CSV must contain one DeltaG column: "
            "'delta_g', 'prodigy_delta_g', 'DeltaG', 'dG', or 'dg'."
        )

    aff = aff[["job_name", dg_col]].copy()
    aff = aff.rename(columns={dg_col: "delta_g"})
    aff["job_name"] = aff["job_name"].astype(str)
    aff["delta_g"] = pd.to_numeric(aff["delta_g"], errors="coerce")
    return aff


def find_structure_for_job(job_dir: str | Path) -> Path | None:
    """
    Find AF3 predicted structure in a job directory.

    Expected possibilities:
      <job_name>_model.cif
      model.cif
      seed-42_sample-0/model.cif
      *.pdb
    """
    if pd.isna(job_dir):
        return None

    p = Path(str(job_dir))

    if p.is_file() and p.suffix.lower() in [".cif", ".pdb"]:
        return p

    if not p.exists() or not p.is_dir():
        return None

    candidates = []

    candidates.extend(sorted(p.glob("*_model.cif")))
    candidates.extend(sorted(p.glob("model.cif")))
    candidates.extend(sorted(p.glob("seed-*_*/*model.cif")))
    candidates.extend(sorted(p.glob("seed-*/model.cif")))
    candidates.extend(sorted(p.glob("*.pdb")))
    candidates.extend(sorted(p.glob("seed-*_*/*.pdb")))
    candidates.extend(sorted(p.glob("seed-*/*.pdb")))

    if candidates:
        return candidates[0]

    return None


def parse_prodigy_delta_g(text: str) -> float:
    """
    Parse PRODIGY binding affinity output.

    Typical lines include variants like:
      Predicted binding affinity (kcal.mol-1):     -10.4
      Predicted binding affinity (kcal/mol):       -10.4
    """
    patterns = [
        r"Predicted binding affinity.*?:\s*([-+]?\d+(?:\.\d+)?)",
        r"binding affinity.*?:\s*([-+]?\d+(?:\.\d+)?)",
        r"DeltaG.*?:\s*([-+]?\d+(?:\.\d+)?)",
        r"ΔG.*?:\s*([-+]?\d+(?:\.\d+)?)",
    ]

    for pat in patterns:
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            return float(m.group(1))

    raise ValueError("Could not parse PRODIGY DeltaG from output.")


def run_prodigy_on_structure(
    structure_path: Path,
    light_chain: str,
    heavy_chain: str,
    antigen_chain: str,
    prodigy_bin: str = "prodigy",
    timeout: int = 120,
) -> float:
    """
    Run PRODIGY on one predicted complex.

    PRODIGY interface chains are given as antibody chains vs antigen chain:
      A,B C
    for light=A, heavy=B, antigen=C.

    The command used is:
      prodigy structure.cif --selection A,B C

    If your local PRODIGY build uses a different CLI syntax, edit cmd below.
    """
    selection_ab = f"{light_chain},{heavy_chain}"
    selection_ag = f"{antigen_chain}"

    cmd = [
        prodigy_bin,
        str(structure_path),
        "--selection",
        selection_ab,
        selection_ag,
    ]

    result = subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )

    output = (result.stdout or "") + "\n" + (result.stderr or "")

    if result.returncode != 0:
        raise RuntimeError(
            f"PRODIGY failed for {structure_path}\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return code: {result.returncode}\n"
            f"Output:\n{output}"
        )

    return parse_prodigy_delta_g(output)


def load_prodigy_cache(path: str | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["job_name", "delta_g", "structure_path", "prodigy_status"])

    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["job_name", "delta_g", "structure_path", "prodigy_status"])

    cache = pd.read_csv(p)
    required = {"job_name", "delta_g"}
    if not required.issubset(set(cache.columns)):
        return pd.DataFrame(columns=["job_name", "delta_g", "structure_path", "prodigy_status"])

    cache["job_name"] = cache["job_name"].astype(str)
    cache["delta_g"] = pd.to_numeric(cache["delta_g"], errors="coerce")
    return cache


def save_prodigy_cache(cache: pd.DataFrame, path: str | None) -> None:
    if path is None:
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    cache.to_csv(path, index=False)


def add_missing_prodigy_delta_g(
    df: pd.DataFrame,
    prodigy_bin: str = "prodigy",
    prodigy_cache: str | None = None,
    structure_col: str = "job_dir",
    light_chain_col: str = "light_chain",
    heavy_chain_col: str = "heavy_chain",
    antigen_chain_col: str = "antigen_chain",
    timeout: int = 120,
) -> pd.DataFrame:
    """
    Fill missing delta_g values by running PRODIGY.

    Requires:
      job_name
      job_dir or structure path column
      light_chain
      heavy_chain
      antigen_chain

    Existing delta_g values are preserved.
    Cached values are used before running PRODIGY.
    """
    df = df.copy()

    if "delta_g" not in df.columns:
        df["delta_g"] = float("nan")

    cache = load_prodigy_cache(prodigy_cache)
    if len(cache) > 0:
        cache_small = cache[["job_name", "delta_g"]].dropna().drop_duplicates("job_name")
        df = df.merge(
            cache_small.rename(columns={"delta_g": "delta_g_cached"}),
            on="job_name",
            how="left",
        )
        missing = df["delta_g"].isna() & df["delta_g_cached"].notna()
        df.loc[missing, "delta_g"] = df.loc[missing, "delta_g_cached"]
        df = df.drop(columns=["delta_g_cached"])

    new_cache_rows = []

    for idx, row in df.iterrows():
        if not pd.isna(row.get("delta_g", float("nan"))):
            continue

        job_name = str(row["job_name"])

        if structure_col not in row or pd.isna(row[structure_col]):
            new_cache_rows.append(
                {
                    "job_name": job_name,
                    "delta_g": float("nan"),
                    "structure_path": "",
                    "prodigy_status": f"missing structure_col {structure_col}",
                }
            )
            continue

        structure_path = find_structure_for_job(row[structure_col])

        if structure_path is None:
            new_cache_rows.append(
                {
                    "job_name": job_name,
                    "delta_g": float("nan"),
                    "structure_path": str(row[structure_col]),
                    "prodigy_status": "structure_not_found",
                }
            )
            continue

        light_chain = str(row.get(light_chain_col, "A"))
        heavy_chain = str(row.get(heavy_chain_col, "B"))
        antigen_chain = str(row.get(antigen_chain_col, "C"))

        try:
            dg = run_prodigy_on_structure(
                structure_path=structure_path,
                light_chain=light_chain,
                heavy_chain=heavy_chain,
                antigen_chain=antigen_chain,
                prodigy_bin=prodigy_bin,
                timeout=timeout,
            )
            df.at[idx, "delta_g"] = dg

            new_cache_rows.append(
                {
                    "job_name": job_name,
                    "delta_g": dg,
                    "structure_path": str(structure_path),
                    "prodigy_status": "ok",
                }
            )

            print(f"[PRODIGY OK] {job_name}: delta_g={dg:.3f}")

        except Exception as e:
            new_cache_rows.append(
                {
                    "job_name": job_name,
                    "delta_g": float("nan"),
                    "structure_path": str(structure_path),
                    "prodigy_status": f"failed: {e}",
                }
            )
            print(f"[PRODIGY FAIL] {job_name}: {e}")

    if new_cache_rows:
        new_cache = pd.DataFrame(new_cache_rows)
        combined = pd.concat([cache, new_cache], ignore_index=True)
        combined = combined.drop_duplicates("job_name", keep="last")
        save_prodigy_cache(combined, prodigy_cache)

    return df


def compute_otr(
    df: pd.DataFrame,
    affinity_df: pd.DataFrame | None = None,
    group_column: str | None = None,
    use_bsa_penalty: bool = False,
) -> pd.DataFrame:
    df = df.copy()

    if "job_name" not in df.columns:
        raise ValueError("Input metrics CSV must contain 'job_name'.")

    df["job_name"] = df["job_name"].astype(str)

    if affinity_df is not None:
        if "delta_g" in df.columns:
            df = df.drop(columns=["delta_g"])
        df = df.merge(affinity_df, on="job_name", how="left")
    elif "delta_g" not in df.columns:
        df["delta_g"] = float("nan")

    if "label" not in df.columns:
        df["label"] = df["job_name"].apply(infer_label)

    if group_column is not None:
        if group_column not in df.columns:
            raise ValueError(f"group_column '{group_column}' not found.")
        df["group_id"] = df[group_column].astype(str)
    elif "group_id" not in df.columns:
        df["group_id"] = df["job_name"].apply(infer_group)

    df["iptm_global_term"] = pd.to_numeric(df.get("iptm"), errors="coerce")
    df["iptm_AbAg_term"] = pd.to_numeric(df.get("iptm_AbAg_mean"), errors="coerce")

    df["pae_AbAg_used"] = df.apply(choose_pae, axis=1)
    df["phi_PAE_AbAg_clean"] = df["pae_AbAg_used"].apply(phi_pae)

    df["delta_g"] = pd.to_numeric(df["delta_g"], errors="coerce")
    df["phi_affinity_clean"] = df["delta_g"].apply(sigmoid_affinity)

    df["P_penalty_clean"] = df.apply(
        lambda row: compute_penalty(row, use_bsa_penalty=use_bsa_penalty),
        axis=1,
    )

    df["OTR_clean"] = (
        0.30 * df["iptm_global_term"]
        + 0.35 * df["iptm_AbAg_term"]
        + 0.20 * df["phi_PAE_AbAg_clean"]
        + 0.15 * df["phi_affinity_clean"]
        - df["P_penalty_clean"]
    )

    df["OTR_clean_noaff"] = (
        0.35 * df["iptm_global_term"]
        + 0.40 * df["iptm_AbAg_term"]
        + 0.25 * df["phi_PAE_AbAg_clean"]
        - df["P_penalty_clean"]
    )

    df["rank_OTR_clean"] = (
        df.groupby("group_id")["OTR_clean"]
        .rank(ascending=False, method="min")
        .astype("Int64")
    )

    df["rank_OTR_clean_noaff"] = (
        df.groupby("group_id")["OTR_clean_noaff"]
        .rank(ascending=False, method="min")
        .astype("Int64")
    )

    margins = []
    for group_id, sub in df.groupby("group_id"):
        pos = sub[sub["label"] == 1]
        dec = sub[sub["label"] == 0]

        pos_score = float("nan")
        best_decoy = float("nan")
        margin = float("nan")

        if len(pos) > 0:
            pos_score = pos["OTR_clean"].max()
        if len(dec) > 0:
            best_decoy = dec["OTR_clean"].max()
        if not pd.isna(pos_score) and not pd.isna(best_decoy):
            margin = pos_score - best_decoy

        margins.append(
            {
                "group_id": group_id,
                "positive_OTR_clean": pos_score,
                "best_decoy_OTR_clean": best_decoy,
                "delta_OTR_clean": margin,
            }
        )

    margin_df = pd.DataFrame(margins)
    df = df.merge(margin_df, on="group_id", how="left")

    return df


def call_selectivity(delta: float) -> str:
    if pd.isna(delta):
        return "not_evaluable"
    if delta > 0.15:
        return "clean_selectivity"
    if delta > 0.05:
        return "moderate_margin"
    if delta >= 0.0:
        return "borderline_off_target_risk"
    return "decoy_scores_above_cognate"


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for group_id, sub in df.groupby("group_id"):
        pos = sub[sub["label"] == 1]
        if len(pos) == 0:
            continue

        pos_row = pos.sort_values("OTR_clean", ascending=False).iloc[0]

        rows.append(
            {
                "group_id": group_id,
                "positive_job": pos_row["job_name"],
                "positive_OTR_clean": pos_row["OTR_clean"],
                "positive_rank_OTR_clean": pos_row["rank_OTR_clean"],
                "delta_OTR_clean": pos_row["delta_OTR_clean"],
                "selectivity_call": call_selectivity(pos_row["delta_OTR_clean"]),
            }
        )

    return pd.DataFrame(rows).sort_values("group_id")


def main():
    parser = argparse.ArgumentParser(
        description="Compute clean OTR score from AF3 engagement metrics and PRODIGY DeltaG."
    )

    parser.add_argument("--metrics_csv", required=True)
    parser.add_argument("--affinity_csv", default=None)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--summary_csv", default=None)
    parser.add_argument("--group_column", default=None)
    parser.add_argument("--use_bsa_penalty", action="store_true")

    parser.add_argument(
        "--run_prodigy_if_missing",
        action="store_true",
        help="Run PRODIGY for rows with missing delta_g.",
    )

    parser.add_argument(
        "--prodigy_bin",
        default="prodigy",
        help="PRODIGY executable name or path.",
    )

    parser.add_argument(
        "--prodigy_cache",
        default="prodigy_delta_g_cache.csv",
        help="CSV cache for PRODIGY DeltaG values.",
    )

    parser.add_argument(
        "--structure_col",
        default="job_dir",
        help="Column containing AF3 job directory or structure path.",
    )

    parser.add_argument(
        "--prodigy_timeout",
        type=int,
        default=120,
        help="Timeout in seconds per PRODIGY job.",
    )

    args = parser.parse_args()

    metrics = pd.read_csv(args.metrics_csv)
    affinity = load_affinity_table(args.affinity_csv)

    if affinity is not None:
        if "delta_g" in metrics.columns:
            metrics = metrics.drop(columns=["delta_g"])
        metrics = metrics.merge(affinity, on="job_name", how="left")

    if args.run_prodigy_if_missing:
        metrics = add_missing_prodigy_delta_g(
            metrics,
            prodigy_bin=args.prodigy_bin,
            prodigy_cache=args.prodigy_cache,
            structure_col=args.structure_col,
            timeout=args.prodigy_timeout,
        )

    scored = compute_otr(
        metrics,
        affinity_df=None,
        group_column=args.group_column,
        use_bsa_penalty=args.use_bsa_penalty,
    )

    scored.to_csv(args.output_csv, index=False)

    if args.summary_csv:
        summary = summarize(scored)
        summary.to_csv(args.summary_csv, index=False)

    print(f"Wrote scored interactions: {args.output_csv}")
    if args.summary_csv:
        print(f"Wrote group summary: {args.summary_csv}")

    print(f"Rows: {len(scored)}")
    print(f"Positives: {int((scored['label'] == 1).sum())}")
    print(f"Decoys: {int((scored['label'] == 0).sum())}")
    print(f"Rows with full OTR_clean: {int(scored['OTR_clean'].notna().sum())}")
    print(f"Rows with OTR_clean_noaff: {int(scored['OTR_clean_noaff'].notna().sum())}")


if __name__ == "__main__":
    main()