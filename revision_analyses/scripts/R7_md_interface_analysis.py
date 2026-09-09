#!/usr/bin/env python3
"""
Revision analysis 7 -- Interface-resolved analysis of the 500 ns scFv-73-0/CLDN4
membrane simulation.  Addresses Reviewer 1 major point 5 and Reviewer 2 major point 7.

The submitted manuscript reported only global backbone RMSD and radius of gyration,
which cannot distinguish a maintained interface from a complex that has drifted while
each partner remains internally folded.  Here we add the interface-resolved metrics
the reviewers requested: interface RMSD, contact number and per-residue contact
occupancy, hydrogen-bond persistence, salt bridges, buried surface area, and
per-chain RMSF.

Chains: PROA = VH (121 aa), PROB = VL (108 aa), PROC = CLDN4 (180 aa, construct
residue 1 = UniProt O14493 residue 5).

Outputs: output/md_interface_timeseries.csv
         output/md_contact_occupancy.csv
         output/md_interface_summary.txt
"""
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array

G = "/mnt/ssd2/claude-tmp/claude-1000/-mnt-ssd1-Projects-DesCARTes-wd/19d8a1ff-bbde-4db4-b38a-efe491578993/scratchpad/gmx/gromacs/"
OFFSET = 4          # construct residue i -> UniProt residue i + OFFSET
CUT = 4.5           # heavy-atom contact cut-off (A)
EQ_NS = 100         # frames before this are treated as equilibration

u = mda.Universe(G + "step7_production_500ns.tpr", G + "step7_pbc_fixed.xtc")
ab = u.select_atoms("segid seg_0_PROA seg_1_PROB")     # scFv (VH + VL)
ag = u.select_atoms("segid seg_2_PROC")                # CLDN4
ab_h = ab.select_atoms("not name H*")
ag_h = ag.select_atoms("not name H*")
times_ns = np.array([ts.time / 1000.0 for ts in u.trajectory])
n = len(u.trajectory)
print(f"{n} frames, {times_ns[0]:.0f}-{times_ns[-1]:.0f} ns, "
      f"{u.atoms.n_atoms} atoms")

# ---- interface residues defined on the first frame ---------------------------
u.trajectory[0]
d0 = distance_array(ab_h.positions, ag_h.positions)
ab_if0 = {ab_h[i].residue.resindex for i in np.where((d0 < CUT).any(axis=1))[0]}
ag_if0 = {ag_h[j].residue.resindex for j in np.where((d0 < CUT).any(axis=0))[0]}
print(f"initial interface: {len(ab_if0)} scFv residues, {len(ag_if0)} CLDN4 residues")
if_sel = u.select_atoms(
    "backbone and resindex " + " ".join(str(i) for i in sorted(ab_if0 | ag_if0)))

# ---- time series -------------------------------------------------------------
ref = u.copy(); ref.trajectory[0]
rows = []
occ_ab, occ_ag, occ_pair = {}, {}, {}
for k, ts in enumerate(u.trajectory):
    d = distance_array(ab_h.positions, ag_h.positions)
    mask = d < CUT
    n_atom_contacts = int(mask.sum())
    ai, aj = np.where(mask)
    pab = {ab_h[i].residue.resindex for i in ai}
    pag = {ag_h[j].residue.resindex for j in aj}
    for r in pab: occ_ab[r] = occ_ab.get(r, 0) + 1
    for r in pag: occ_ag[r] = occ_ag.get(r, 0) + 1
    for i, j in zip(ai, aj):
        key = (ab_h[i].residue.resindex, ag_h[j].residue.resindex)
        occ_pair[key] = occ_pair.get(key, 0) + 1
    rows.append(dict(time_ns=ts.time / 1000.0,
                     atom_contacts=n_atom_contacts,
                     residue_contacts_scFv=len(pab),
                     residue_contacts_CLDN4=len(pag),
                     residue_pairs=len({(ab_h[i].residue.resindex,
                                         ag_h[j].residue.resindex)
                                        for i, j in zip(ai, aj)})))
ts_df = pd.DataFrame(rows)

# ---- RMSDs (antigen-superposed = interface-relevant) -------------------------
R = rms.RMSD(u, ref,
             select="segid seg_2_PROC and backbone",                 # superpose on CLDN4
             groupselections=["segid seg_0_PROA seg_1_PROB and backbone",  # scFv drift
                              "protein and backbone",
                              "backbone and resindex " +
                              " ".join(str(i) for i in sorted(ab_if0 | ag_if0))])
R.run()
ts_df["CLDN4_backbone_rmsd"] = R.results.rmsd[:, 2]
ts_df["scFv_rmsd_antigen_aligned"] = R.results.rmsd[:, 3]
ts_df["complex_backbone_rmsd"] = R.results.rmsd[:, 4]
ts_df["interface_rmsd"] = R.results.rmsd[:, 5]

# ---- hydrogen bonds across the interface ------------------------------------
hb = HydrogenBondAnalysis(
    universe=u,
    donors_sel=None,
    between=["segid seg_0_PROA seg_1_PROB", "segid seg_2_PROC"],
    d_a_cutoff=3.5, d_h_a_angle_cutoff=150, update_selections=False)
hb.run()
counts = hb.count_by_time()
ts_df["interface_hbonds"] = counts
ts_df.to_csv("output/md_interface_timeseries.csv", index=False)

# ---- salt bridges ------------------------------------------------------------
pos = u.select_atoms("(resname ARG LYS and name NH* NZ NE) or "
                     "(resname HSP HIS and name ND1 NE2)")
neg = u.select_atoms("(resname ASP GLU and name OD* OE*)")
sb_count = []
for ts in u.trajectory:
    a = pos.select_atoms("segid seg_0_PROA seg_1_PROB")
    b = neg.select_atoms("segid seg_2_PROC")
    c = pos.select_atoms("segid seg_2_PROC")
    e = neg.select_atoms("segid seg_0_PROA seg_1_PROB")
    k = 0
    if len(a) and len(b): k += int((distance_array(a.positions, b.positions) < 4.0).any(axis=1).sum())
    if len(c) and len(e): k += int((distance_array(c.positions, e.positions) < 4.0).any(axis=1).sum())
    sb_count.append(k)
ts_df["salt_bridges"] = sb_count

# ---- buried surface area (every 5th frame) -----------------------------------
try:
    import freesasa
    def sasa(atoms):
        rad = [1.9 if a.name.startswith("C") else 1.7 if a.name.startswith("N")
               else 1.52 if a.name.startswith("O") else 1.8 for a in atoms]
        s = freesasa.calcCoord(atoms.positions.flatten().tolist(), rad)
        return s.totalArea()
    bsa = []
    for ts in u.trajectory[::5]:
        bsa.append(dict(time_ns=ts.time / 1000.0,
                        bsa=(sasa(ab_h) + sasa(ag_h) - sasa(ab_h + ag_h)) / 2.0))
    pd.DataFrame(bsa).to_csv("output/md_bsa.csv", index=False)
    bsa_v = np.array([b["bsa"] for b in bsa]); bsa_t = np.array([b["time_ns"] for b in bsa])
except Exception as exc:
    print("freesasa unavailable:", exc); bsa_v = None

# ---- per-residue occupancy ---------------------------------------------------
def rec(occ, tag):
    out = []
    for ridx, c in occ.items():
        r = u.residues[ridx]
        num = r.resid + OFFSET if tag == "CLDN4" else r.resid
        out.append(dict(chain=tag, resid_construct=r.resid, resid_uniprot=num,
                        resname=r.resname, occupancy=c / n))
    return out
occ_df = pd.DataFrame(rec(occ_ab, "scFv") + rec(occ_ag, "CLDN4"))
occ_df.sort_values(["chain", "occupancy"], ascending=[True, False]) \
      .round(3).to_csv("output/md_contact_occupancy.csv", index=False)
ts_df.round(3).to_csv("output/md_interface_timeseries.csv", index=False)

# ---- report ------------------------------------------------------------------
L = []
def emit(s=""):
    print(s); L.append(s)

post = ts_df[ts_df.time_ns >= EQ_NS]
emit("=" * 84)
emit("INTERFACE-RESOLVED ANALYSIS OF THE 500 ns scFv-73-0/CLDN4 SIMULATION")
emit("=" * 84)
emit(f"Frames: {n} (1 ns spacing). Summary window: {EQ_NS}-500 ns (n = {len(post)}).")
emit("")
def stat(col, unit=""):
    v = post[col]
    emit(f"  {col:32s} mean {v.mean():8.2f} +/- {v.std():5.2f} {unit}"
         f"   range {v.min():.2f}-{v.max():.2f}")
emit("Structural stability:")
for c in ["complex_backbone_rmsd", "CLDN4_backbone_rmsd",
          "scFv_rmsd_antigen_aligned", "interface_rmsd"]:
    stat(c, "A")
emit("")
emit("Interface persistence:")
for c in ["atom_contacts", "residue_pairs", "interface_hbonds", "salt_bridges"]:
    stat(c)
if bsa_v is not None:
    m = bsa_t >= EQ_NS
    emit(f"  {'buried surface area':32s} mean {bsa_v[m].mean():8.1f} +/- "
         f"{bsa_v[m].std():5.1f} A^2   range {bsa_v[m].min():.0f}-{bsa_v[m].max():.0f}")
emit("")
frac_lost = float((post.atom_contacts == 0).mean())
emit(f"Fraction of the summary window with zero interface contacts: {frac_lost:.3f}")
emit(f"Contacts at 500 ns relative to the mean over 100-500 ns: "
     f"{ts_df.atom_contacts.iloc[-1] / post.atom_contacts.mean():.2f}")
emit("")
emit("Most persistent CLDN4 epitope residues (occupancy = fraction of 501 frames):")
top = occ_df[occ_df.chain == "CLDN4"].sort_values("occupancy", ascending=False).head(15)
for _, r in top.iterrows():
    emit(f"  {r.resname}{int(r.resid_uniprot):<4d} (construct {int(r.resid_construct):<4d})"
         f"  occupancy {r.occupancy:.2f}")
emit("")
emit("Most persistent scFv paratope residues:")
top2 = occ_df[occ_df.chain == "scFv"].sort_values("occupancy", ascending=False).head(12)
for _, r in top2.iterrows():
    emit(f"  {r.resname}{int(r.resid_construct):<4d}  occupancy {r.occupancy:.2f}")

open("output/md_interface_summary.txt", "w").write("\n".join(L) + "\n")
print("\nSaved -> output/md_interface_{timeseries,summary}, md_contact_occupancy.csv")
