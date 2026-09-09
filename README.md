# DesCARTes — STRIVE

Code and model outputs for:

> **STRIVE: An Integrated In Silico Framework for Structure-Guided De Novo Design and
> Off-Target Risk Prioritization of CLDN4-Targeting scFv Candidates**

Archived release: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22680234.svg)](https://doi.org/10.5281/zenodo.22680234)

STRIVE (Structural Therapeutic Risk and Interaction Viability Engine) is a **model-derived
prioritization score**, not an affinity or safety predictor. Its measured discrimination on
an independent panel of solved antibody–antigen complexes is weak (AUROC 0.606, 95% CI
0.348–0.846) and it does not outperform interface ipTM used alone. Please read
`revision_analyses/output/benchmark_discrimination.csv` before reusing it.

## Layout

```
pipeline/
  design/         RFantibody (RFdiffusion + ProteinMPNN + RF2) driver scripts, PRODIGY, AF3 input prep
  alphafold3/     Nextflow workflow for AF3 MSA + inference; AF3 output parser
  scoring/        STRIVE implementation and the refined/clean OTR variants
  md/             GROMACS .mdp inputs, topology, run scripts, and analysis .xvg outputs
expression/       CLDN4 target selection: TCGA/LuCA pseudobulk, CPS scoring, HPA and IHC
notebooks/        Figure-generating notebooks
results/
  scoring_tables/ on-target and off-target STRIVE scores, RF2 scores for all 1,000 designs,
                  PRODIGY ΔG tables, clinical monoclonal panel
  benchmark/      independent decoy-panel AF3 metrics and scored outputs
revision_analyses/
  scripts/        R1–R8: every number added during peer review, reproducible end to end
  output/         their outputs (specificity margins, benchmark, seed variability,
                  paired IHC statistics, LuCA FDR statistics, HPA pan-tissue, MD interface)
figures/          Supplementary Figure 3 (vector)
environment/      pinned Python environment for revision_analyses
SOFTWARE_VERSIONS.md   full provenance: versions, model checkpoints, database snapshots
MANIFEST.md            what each file is
```

## Reproducing the revision analyses

```bash
pip install -r environment/analysis_requirements.txt
cd revision_analyses
python3 scripts/R1_specificity_reanalysis.py     # Table 1, worst-case margins
python3 scripts/R2_strive_benchmark.py           # decoy benchmark + ablation
python3 scripts/R3_af3_seed_variability.py       # AF3 seed-to-seed variability
python3 scripts/R4_ihc_paired_stats.py           # paired IHC statistics
python3 scripts/R5_luca_pseudobulk_stats.py      # LuCA extraction (needs the atlas)
python3 scripts/R5b_luca_stats_from_profiles.py  # FDR + donor-level statistics
python3 scripts/R6_hpa_pan_normal_tissue.py      # HPA pan-normal-tissue / pan-cancer
python3 scripts/R7_md_interface_analysis.py      # MD interface metrics (needs trajectory)
python3 scripts/R8_supp_fig3_mab_heatmap.py      # Supplementary Figure 3
```

R5 needs the LuCA atlas (`202603_luca.h5ad`) and R7 needs the 500 ns trajectory; both are
too large for git and are archived in the Zenodo record.

## Not in this repository

Kept in the Zenodo record because of size: the 500 ns production trajectory
(`step7_production_500ns.xtc`), the AlphaFold3 model outputs for all 975 antibody–antigen
pairs, the CHARMM-GUI `toppar/` force-field bundle, and the LuCA single-cell atlas.

## Reproducibility caveats

Stated plainly, because they affect what this code can be used for:

1. The AlphaFold3 container was referenced as `alphafold3:latest`, not by digest, so the
   exact build used for the April 2026 screens cannot be recovered. See `SOFTWARE_VERSIONS.md`.
2. Each antibody–antigen pair was predicted with a **single seed and a single diffusion
   sample**. Measured run-to-run variability is ~0.03 STRIVE units (max 0.068); differences
   smaller than ~0.07 units should not be interpreted.
3. The 1,000 designed scFv sequences are not yet deposited here as a single FASTA — see
   `MANIFEST.md`.

## Citation

Please cite the paper and the archived release (DOI 10.5281/zenodo.22680234).
