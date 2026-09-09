# Manifest

Generated for the archived release. Sizes in KB.

> **Note on the 1,000 designed sequences.** Reviewer 2 asked for all 1,000 designed scFv
> sequences to be deposited. `results/scoring_tables/3_rf2.sc` contains the RoseTTAFold2
> metrics for all 1,000 designs keyed by design tag, but the sequences themselves live in
> the RFantibody design PDB outputs and are not yet extracted into a single FASTA. This is
> an outstanding deposit item.


## `pipeline/design/`

RFantibody driver scripts (RFdiffusion backbone generation, ProteinMPNN sequence design, RF2 filtering), PRODIGY runner, AF3 input preparation, clinical mAb combinatorial builder.

- `alternative_rfantibody_slurm.sh` — 6.8 KB
- `antibody_full_pipeline.sh` — 5.0 KB
- `mab_offtarget_af3_combinator.py` — 12.6 KB
- `pdb_to_af3.py` — 5.4 KB
- `run_prodigy.sh` — 0.5 KB

## `pipeline/alphafold3/`

Nextflow workflow splitting AF3 MSA generation from inference, plus the AF3 confidence-output parser used to build STRIVE inputs.

- `af3_otr_clean.nf` — 9.2 KB
- `af3_otr_pipeline.py` — 35.3 KB
- `main.nf` — 2.6 KB
- `nextflow.config` — 1.2 KB
- `run_mab_nf_otrtest.sh` — 0.2 KB
- `run_mab_nf_otrtest_99.sh` — 0.2 KB
- `run_mab_nf_otrtest_99_clean.sh` — 0.2 KB

## `pipeline/scoring/`

STRIVE scoring implementations (original, 'clean' and literature-refined variants) and CDR range definitions.

- `cdr_ranges_default_refined.csv` — 0.1 KB
- `compute_otr_clean.py` — 17.2 KB
- `compute_otr_refined_lit.py` — 15.4 KB
- `subroutines.py` — 1.2 KB

## `pipeline/md/`

GROMACS inputs for the 500 ns membrane simulation: six-step CHARMM-GUI equilibration .mdp files, production .mdp, topology, run scripts, and the RMSD/Rg/RMSF .xvg outputs.

- `MD_equilibration.sh` — 2.7 KB
- `MD_production_500ns.sh` — 3.5 KB
- `MD_trajectory.sh` — 1.2 KB
- `mdout.mdp` — 10.9 KB
- `rg_fixed.xvg` — 24.4 KB
- `rg_scfv73_CLDN4.xvg` — 24.4 KB
- `rg_scfv73_CLDN4_clean.xvg` — 24.4 KB
- `rmsd_backbone_fixed.xvg` — 13.5 KB
- `rmsd_scfv73_CLDN4.xvg` — 13.5 KB
- `rmsd_scfv73_CLDN4_clean.xvg` — 14.3 KB
- `rmsf_chain_A.xvg` — 2.5 KB
- `rmsf_chain_B.xvg` — 2.3 KB
- `rmsf_chain_C.xvg` — 3.4 KB
- `rmsf_fixed.xvg` — 6.7 KB
- `step6.0_minimization.mdp` — 0.6 KB
- `step6.1_equilibration.mdp` — 1.1 KB
- `step6.2_equilibration.mdp` — 1.1 KB
- `step6.3_equilibration.mdp` — 1.3 KB
- `step6.4_equilibration.mdp` — 1.3 KB
- `step6.5_equilibration.mdp` — 1.3 KB
- `step6.6_equilibration.mdp` — 1.3 KB
- `step7_production.mdp` — 1.1 KB
- `topol.top` — 1.0 KB

## `expression/`

CLDN4 target characterisation: TCGA and LuCA pseudobulk analyses, CPS target-selection scoring, HPA immunohistochemistry extraction, and the matched-cohort IHC analysis.

- `20260313_cldn4_IHC_our_data.py` — 4.3 KB
- `20260313_cldn4_hpa.R` — 8.1 KB
- `20260313_cldn4_v24_hpa.R` — 8.9 KB
- `CPS.v1.py` — 11.9 KB
- `CPS.v2.py` — 7.8 KB
- `CPS.v3.py` — 17.2 KB
- `Combining_with_surfaceome.Rmd` — 1.2 KB
- `Downsample-Pseudobulk-DEG.py` — 13.1 KB
- `Luca_Epithelial_DEG.Rmd` — 4.2 KB
- `cldn4_hpa.R` — 8.1 KB
- `cldn4_ihc_ourown.py` — 4.3 KB
- `compute.Brain.expr.py` — 1.9 KB

## `notebooks/`

Figure-generating Jupyter notebooks for the STRIVE screen, the clinical mAb panel, MD RMSD/Rg, and HPA/IHC visualisation.

- `20260310_luca_cldn4_ebru.ipynb` — 5457.9 KB
- `HPA_IHC_visual.ipynb` — 102.0 KB
- `OTR_offtarget.ipynb` — 409.7 KB
- `RFantibody_scoring_visualization.ipynb` — 1744.2 KB
- `RMSD_RG.ipynb` — 230.4 KB
- `decision.ipynb` — 499.4 KB
- `heatmap.ipynb` — 852.1 KB
- `mab_OTR.ipynb` — 1329.7 KB
- `scfv_specificity_decision_notebook.ipynb` — 699.4 KB

## `results/scoring_tables/`

STRIVE scores: 39 candidates vs CLDN4 (on-target) and vs 24 homologues (off-target), RF2 metrics for all 1,000 designs, PRODIGY dG tables, and the clinical monoclonal panel.

- `3_rf2.sc` — 98.5 KB
- `OTR_offtarget_test1.csv` — 77.7 KB
- `OTR_ontarget.csv` — 4.2 KB
- `mab_offtarget_scores.csv` — 68.5 KB
- `mab_prodigy.tsv` — 10.2 KB
- `off_target_prodigy.tsv` — 42.6 KB
- `ontarget_af3_prodigy.tsv` — 1.5 KB

## `results/benchmark/`

Independent decoy-panel benchmark: AF3 metrics and scored outputs for 9 cognate antibody-antigen pairs and 90 non-cognate decoys.

- `af3_otr_metrics.README.txt` — 0.6 KB
- `af3_otr_metrics_99.README.txt` — 0.6 KB
- `af3_otr_metrics_99.csv` — 98.0 KB
- `af3_otr_metrics_99_freesasa_cdr.README.txt` — 0.6 KB
- `af3_otr_metrics_99_freesasa_proxy.README.txt` — 0.6 KB
- `otr_clean_scored_99.csv` — 120.3 KB
- `otr_clean_summary_99.csv` — 1.0 KB
- `otr_refined_lit_99.csv` — 143.4 KB
- `otr_refined_lit_summary_99.csv` — 1.5 KB
- `prodigy_delta_g_cache_99.csv` — 19.1 KB

## `revision_analyses/scripts/`

Analyses added during peer review; each reproduces specific reported values.

- `R1_specificity_reanalysis.py` — 5.6 KB
- `R2_strive_benchmark.py` — 4.8 KB
- `R3_af3_seed_variability.py` — 3.5 KB
- `R4_ihc_paired_stats.py` — 4.1 KB
- `R5_luca_pseudobulk_stats.py` — 6.1 KB
- `R5b_luca_stats_from_profiles.py` — 4.3 KB
- `R6_hpa_pan_normal_tissue.py` — 5.4 KB
- `R7_md_interface_analysis.py` — 8.8 KB
- `R8_supp_fig3_mab_heatmap.py` — 7.1 KB

## `revision_analyses/output/`

Outputs of the above.

- `STRIVE_matrix_39x25.csv` — 18.5 KB
- `Table1_revised.csv` — 2.1 KB
- `TableS_CLDN4_normal_tissue.csv` — 7.8 KB
- `TableS_CLDN4_pan_cancer.csv` — 1.2 KB
- `TableS_IHC_pairs.csv` — 0.7 KB
- `TableS_specificity_full.csv` — 5.2 KB
- `af3_seed_variability.csv` — 0.6 KB
- `af3_seed_variability_long.csv` — 0.8 KB
- `benchmark_discrimination.csv` — 0.4 KB
- `benchmark_per_group.csv` — 0.7 KB
- `benchmark_per_pair.csv` — 125.6 KB
- `cldn4_hpa_normal_all.csv` — 10.3 KB
- `cldn4_hpa_pathology_all.csv` — 1.5 KB
- `hpa_normal_tissue_summary.txt` — 6.6 KB
- `ihc_paired_stats.txt` — 1.0 KB
- `luca_cldn4_profiles.csv` — 188.1 KB
- `luca_report.txt` — 4.8 KB
- `luca_report_FINAL.txt` — 2.6 KB
- `luca_stage_stats.csv` — 1.3 KB
- `luca_stage_stats_FINAL.csv` — 1.6 KB
- `mab_panel_cognate_ranks.csv` — 0.4 KB
- `mab_panel_metadata.csv` — 4.8 KB
- `md_bsa.csv` — 2.4 KB
- `md_contact_occupancy.csv` — 2.2 KB
- `md_contact_occupancy_CLDN4_annotated.csv` — 1.5 KB
- `md_interface_summary.txt` — 3.0 KB
- `md_interface_timeseries.csv` — 23.1 KB
- `worstcase_margins.csv` — 36.3 KB

## `figures/`

Supplementary Figure 3, regenerated as vector artwork.

- `SupplementaryFigure3_mAb_panel.pdf` — 28.3 KB

## `environment/`

Pinned Python environment for the revision analyses.

- `analysis_requirements.txt` — 0.3 KB
