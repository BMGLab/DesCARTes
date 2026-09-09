# Software versions, model checkpoints and database snapshots

Recorded from the analysis hosts on 9 September 2026. Everything here is also stated in
Methods §2.6 of the manuscript.

## Design (de novo scFv generation and filtering)

| Component | Version / identifier | How it was determined |
|---|---|---|
| RFantibody pipeline | `RosettaCommons/RFantibody` commit `184ea07c558797a5c40838b1a9ee7ea3f940f51c` (2 Oct 2025) | `git log -1` in the local clone |
| RFdiffusion checkpoint | `RFdiffusion_Ab.pt` (461 MB) | distributed with RFantibody |
| ProteinMPNN checkpoint | `ProteinMPNN_v48_noise_0.2.pt` (6.4 MB) | distributed with RFantibody |
| RoseTTAFold2 (antibody) checkpoint | `RF2_ab.pt` (281 MB) | distributed with RFantibody |
| Python (pipeline runtime) | 3.10 | RFantibody virtualenv |
| PyTorch | 2.3.1+cu118 | virtualenv dist-info |
| DGL | 2.4.0+cu118 | virtualenv dist-info |
| e3nn | 0.5.1 | virtualenv dist-info |
| NumPy / SciPy | 1.26.4 / 1.13.1 | virtualenv dist-info |
| PyMOL | 3.0 | antigen preparation |

ProteinMPNN sampling temperature 0.2, 10 sequences per backbone, 100 backbones
(1,000 designs). CDR loop-length ranges L1 8–13, L2 7, L3 9–11, H1 7, H2 6, H3 5–13.
Framework held fixed on hu-4D5-8-Fv (HLT format). Filters: target-aligned antibody
RMSD ≤ 5 Å and PAE ≤ 10.

## Structure prediction and scoring

| Component | Version / identifier | How it was determined |
|---|---|---|
| AlphaFold3 | Docker image `alphafold3:latest` — **tag not pinned** | `nextflow.config` |
| AlphaFold3 (host image at time of writing) | v3.0.2, source commit `97639fff6fb22c0d9765089026fe296ee506b60a`, image built 2026-05-18 | `version.py` inside the container |
| AlphaFold3 input dialect | `alphafold3`, version 2 | `*_data.json` |
| Random seed | `modelSeeds = [42]`, one diffusion sample per pair | `*_data.json` |
| Nextflow | 24.10.5, build 5935 | `.nextflow.log` |
| PRODIGY | `prodigy-prot` 2.2.5 | pip metadata |

### Known reproducibility limitation

The AlphaFold3 container was referenced by the mutable tag `alphafold3:latest` rather than
by an immutable digest. The screens were run 23–28 April 2026; the image now present on the
host was rebuilt on 2026-05-18. **The exact AlphaFold3 build used at run time therefore
cannot be recovered retrospectively.** Future runs should pin the image digest. This is
stated in the manuscript rather than papered over.

### Sequence database snapshots (AlphaFold3 data pipeline)

Read from `/mnt/ssd0/alphafold3_db`; the dates are encoded in the filenames.

| Database | Snapshot |
|---|---|
| UniRef90 | `uniref90_2022_05.fa` (2022-05) |
| MGnify clusters | `mgy_clusters_2022_05.fa` (2022-05) |
| BFD | `bfd-first_non_consensus_sequences.fasta` |
| UniProt | `uniprot_all_2021_04.fa` (2021-04) |
| PDB seqres (templates) | `pdb_seqres_2022_09_28.fasta` (2022-09-28) |
| PDB mmCIF template set | `mmcif_files/` (2024-10-11) |

Note: the originally submitted Methods stated that MSAs were built with "Jackhmmer and
HHblits against UniRef90, BFD, MGnify and PDB70". That was incorrect. The AlphaFold3 data
pipeline uses jackhmmer for protein MSAs and hmmbuild/hmmsearch against `pdb_seqres` for
templates; it does not use HHblits and does not use PDB70. The Methods have been corrected.

## Molecular dynamics

| Component | Version | How it was determined |
|---|---|---|
| GROMACS | 2024.3 | `mdrun` log header |
| System builder | CHARMM-GUI v3.7 (Membrane Builder) | `gromacs/README` |
| Force field / water | CHARMM36m / TIP3P | `topol.top` |
| MDAnalysis | 2.10.0 | analysis host |
| FreeSASA | 2.2.1 | analysis host |

Production: 250,000,000 steps × 2 fs = **500 ns**, frames every 1 ns (501 frames),
284.4 ns/day. Thermostat v-rescale, τ_t = 1.0 ps, 310.15 K. Barostat C-rescale,
**semi-isotropic**, τ_p = 5.0 ps, 1.0 bar, compressibility 4.5 × 10⁻⁵ bar⁻¹.
System 80,717 atoms; 203 lipids; 17,149 TIP3P waters; 61 K⁺ / 45 Cl⁻ (≈0.15 M KCl).

## Expression and statistical analyses

| Component | Version |
|---|---|
| Python | 3.13.9 (revision analyses); 3.10.12 / 3.12.12 (original figure notebooks) |
| NumPy / Pandas / SciPy | 2.3.5 / 2.3.3 / 1.16.3 |
| statsmodels | 0.14.5 |
| Scanpy / AnnData | 1.11.5 / 0.12.6 |
| Matplotlib / Seaborn | 3.10.7 / 0.13.2 |
| openpyxl | 3.1.5 |
| HPAanalyze (R) | 1.24.0, bundling HPA v22.0 |

## Hardware

- **AlphaFold3**: NVIDIA RTX A4000 (16 GB), 52 CPU threads.
- **Molecular dynamics**: NVIDIA A100 (80 GB), Intel Xeon Platinum 8362 (64 cores).
