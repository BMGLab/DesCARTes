import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gc
import ctypes
import scipy.sparse
from joblib import Parallel, delayed
from scipy.sparse import issparse
import multiprocessing
multiprocessing.set_start_method("spawn", force=True)
from collections import defaultdict


#In case the adata of pseudo-bulk data needs to be saved.

def savePBdata2(adata, filename):
    # Convert counts to dense matrix if it's sparse
    if scipy.sparse.issparse(adata.X):
        counts = adata.X.toarray()
    else:
        counts = adata.X

    # Convert to DataFrame
    counts_df = pd.DataFrame(counts, index=adata.obs_names, columns=adata.var_names)

    # (Optional) Convert to integer
    counts_df = counts_df.astype(int, errors='ignore')

    # Save AnnData object
    adata.write(filename)
    
    # Save counts as CSV for inspection (optional)
    counts_df.to_csv(filename.replace(".h5ad", "_counts.csv"))

    # print(f"Saved AnnData to {filename} and count matrix to CSV.")

import os
import pickle

def parse_gtf(gtf_path):
    """
    Parse a GTF file and return a list of (gene_id, gene_name) tuples.
    Caches the result to a .pkl file to avoid reparsing on future runs.

    Parameters:
        gtf_path (str): Path to the GTF file.

    Returns:
        list of tuples: Each tuple contains (gene_id, gene_name)
    """
    cache_path = gtf_path + ".parsed.pkl"

    # Load from cache if available
    if os.path.exists(cache_path):
        with open(cache_path, 'rb') as f:
            # print(f"Loading cached GTF data from {cache_path}")
            return pickle.load(f)

    # Else parse the GTF
    # print(f"Parsing GTF file {gtf_path}")
    with open(gtf_path, 'r') as f:
        gtf_lines = list(f)

    gtf_filtered = [
        line for line in gtf_lines
        if not line.startswith('#') and 'gene_id "' in line and 'gene_name "' in line
    ]

    gtf_tuples = [
        (
            line.split('gene_id "')[1].split('"')[0],
            line.split('gene_name "')[1].split('"')[0]
        )
        for line in gtf_filtered
    ]

    # Save to cache
    with open(cache_path, 'wb') as f:
        pickle.dump(gtf_tuples, f)
        # print(f"Cached parsed GTF data to {cache_path}")

    return gtf_tuples


def bulk_one_group(adata, group_key, idx, metacols, cell_type_colname, min_cells, max_cells, rng):
    s, ct = group_key
    cell_indices = np.array(idx)
    n_cells_total = len(cell_indices)

    if n_cells_total < min_cells:
        return None

    selected = rng.choice(cell_indices, size=min(max_cells, n_cells_total), replace=False)
    sub = adata[selected]

    summed_X = sub.X.sum(axis=0)
    if issparse(sub.X):
        summed_X = summed_X.A1

    adata_rep = sc.AnnData(X=summed_X.reshape(1, -1), var=sub.var[[]])
    adata_rep = adata_rep.copy()
    adata_rep.obs_names = [f"{s}_{ct}"]
    adata_rep.obs[cell_type_colname] = ct

    for col in metacols:
        adata_rep.obs[col] = sub.obs[col].iloc[0]

    adata_rep.obs["n_cells"] = len(selected)
    adata_rep.obs["sample"] = s
    adata_rep.obs["log2_n_cells"] = np.log2(len(selected) + 1)  # Adding 1 to avoid log(0)
    adata_rep.obs["Contrast"] = ct
    # print(f"Sample {s}, cell type {ct} has {len(selected)} cells, pseudobulk shape: {adata_rep.shape}")
    if adata_rep.shape[1] == 0:
        # print(f"Warning: No genes in pseudobulk for sample {s}, cell type {ct}. Skipping.")
        return None
    
    return adata_rep


def make_pseudobulk_parallel(adata, sample_colname, cell_type_colname, metacols, min_cells=40, max_cells=50, random_state=0, n_jobs=4):
    rng = np.random.default_rng(random_state)
    groups = adata.obs.groupby([sample_colname, cell_type_colname])

    results = Parallel(n_jobs=n_jobs, backend="threading")(
        delayed(bulk_one_group)(
            adata, (s, ct), idx, metacols, cell_type_colname, min_cells, max_cells,
            np.random.default_rng(random_state + i)
        )
        for i, ((s, ct), idx) in enumerate(groups.indices.items())
    )

    ps_list = [res for res in results if res is not None]
    ps_adata = sc.concat(ps_list)
    ps_adata = ps_adata.copy()
    ps_adata.obs["n_cells"] = ps_adata.obs["n_cells"].astype(int)
    ps_adata.obs["log2_n_cells"] = np.log2(ps_adata.obs["n_cells"] + 1)
    savePBdata2(ps_adata, 'results/pseudobulk_data.h5ad')

    return ps_adata


def run_deseq_for_celltype(ps_adata, c, cell_type_colname, cell_of_interest, metacols, gtf, condition_col=['Contrast'], latent_factors=['assay']):
    sub_c = ps_adata[ps_adata.obs[cell_type_colname] == c ]
    sub_c = sub_c.copy()
    sub_c.obs[condition_col] = c
    sub_ct = ps_adata[ps_adata.obs[cell_type_colname] == cell_of_interest ]
    sub_ct = sub_ct.copy()
    sub_ct.obs[condition_col] = cell_of_interest
    pb = sc.concat([sub_c, sub_ct])

    counts = pd.DataFrame(pb.X.toarray() if hasattr(pb.X, 'toarray') else pb.X, columns=pb.var_names)
    counts = counts.astype(int, errors='ignore')
    counts += 1  # Avoid 0s
    
    # Create a design factors:
    if "log2_n_cells" not in latent_factors:
        latent_factors.append("log2_n_cells")
    
    if condition_col not in latent_factors:
        design_factors_local = latent_factors + condition_col
    else:
        design_factors_local = latent_factors[:]
    
    # print("design_factors_local:", f"~ {' + '.join(design_factors_local)}")

    dds = DeseqDataSet(counts=counts, 
                       metadata=pb.obs, 
                       design=f"~ {' + '.join(design_factors_local)}", 
                       quiet=True,
                       low_memory=True)
    sc.pp.filter_genes(dds, min_cells=10)
    # print(f"Running DESeq2 for cell type {c} against {cell_of_interest} with design factors: {design_factors_local}")
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=('Contrast', cell_of_interest, c))
    stat_res.summary()

    de = stat_res.results_df
    de['gene_symbols'] = de.index
    de['gene_name'] = de['gene_symbols'].map(dict(gtf)) if gtf else de['gene_symbols']
    de = de[(de['padj'] < 0.05) & (de['log2FoldChange'] > 1.0)]
    de['SubType'] = c
    de['inContrastTo'] = 'others'
    return de

def run_deseq_for_celltype_v2(
    ps_adata,
    c,
    cell_type_colname,
    cell_of_interest,
    metacols,
    gtf,
    condition_col=['Contrast'],
    latent_factors=['assay']
):
    sub_c = ps_adata[ps_adata.obs[cell_type_colname] == c]
    sub_ct = ps_adata[ps_adata.obs[cell_type_colname] == cell_of_interest]

    # Check that both are non-empty
    if sub_c.n_obs == 0 or sub_ct.n_obs == 0:
        print(f"⚠️  Skipping {c} vs {cell_of_interest} — one group has no cells.")
        return None

    sub_c = sub_c.copy()
    sub_ct = sub_ct.copy()

    sub_c.obs[condition_col] = c
    sub_ct.obs[condition_col] = cell_of_interest

    pb = sc.concat([sub_c, sub_ct])

    if pb.n_obs == 0:
        print(f"⚠️  Empty pseudobulk for {c} vs {cell_of_interest}. Skipping.")
        return None

    counts = (
        pd.DataFrame(pb.X.toarray(), columns=pb.var_names)
        if hasattr(pb.X, 'toarray') else
        pd.DataFrame(pb.X, columns=pb.var_names)
    )
    counts = counts.astype(int, errors='ignore')
    counts += 1  # Pseudo-count to avoid zeros

    if "log2_n_cells" not in latent_factors:
        latent_factors.append("log2_n_cells")

    if condition_col not in latent_factors:
        design_factors_local = latent_factors + condition_col
    else:
        design_factors_local = latent_factors[:]

    dds = DeseqDataSet(
        counts=counts,
        metadata=pb.obs,
        design=f"~ {' + '.join(design_factors_local)}",
        quiet=True,
        low_memory=True,
        n_cpus=1        # default = None
    )

    sc.pp.filter_genes(dds, min_cells=10)

    dds.deseq2()

    # ✅ Check that contrast levels exist:
    design_levels = pb.obs['Contrast'].unique()
    if cell_of_interest not in design_levels or c not in design_levels:
        print(f"⚠️  Contrast levels missing in design: {design_levels}. Skipping {c} vs {cell_of_interest}")
        return None

    stat_res = DeseqStats(dds, contrast=('Contrast', cell_of_interest, c))
    stat_res.summary()

    de = stat_res.results_df
    de['gene_symbols'] = de.index
    de['gene_name'] = de['gene_symbols'].map(dict(gtf)) if gtf else de['gene_symbols']
    #de = de[(de['padj'] < 0.05) & (de['log2FoldChange'] > 1.0)] # report all DEGs
    de['SubType'] = c
    de['inContrastTo'] = cell_of_interest
    return de

def repeat_deseq_pseudobulk_analysis_parallel(adata, ctyps, sample_colname, cell_type_colname, cell_of_interest, metacols, outCsv, N=10, min_cells=40, max_cells=50, min_repeats=5, gtf=None, n_jobs=4):
    gene_hits = defaultdict(int)
    final_results = []

    for i in range(N):
        # print(f"=== Iteration {i+1}/{N} ===")
        ps_adata = make_pseudobulk_parallel(
            adata, 
            sample_colname=sample_colname,
            cell_type_colname=cell_type_colname,
            metacols=metacols,
            min_cells=min_cells,
            max_cells=max_cells,
            random_state=i,
            n_jobs=n_jobs
        )

        results = Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(run_deseq_for_celltype)(ps_adata, c, cell_type_colname, cell_of_interest, metacols, gtf) 
            for c in ctyps
        )
        # print(f"Results from iteration {i+1}: {len(results)} cell types processed.")
        for de in results:
            if de is not None and not de.empty:
                for g, c in zip(de['gene_name'], de['SubType']):
                    gene_hits[(g, c)] += 1
                final_results.append(de)
        # print(f"Total unique genes across all iterations: {len(gene_hits)}")
    all_deg_df = pd.concat(final_results)
    stable_genes = [(gene, subtype) for (gene, subtype), count in gene_hits.items() if count >= min_repeats]
    stable_df = all_deg_df[all_deg_df.apply(lambda row: (row['gene_name'], row['SubType']) in stable_genes, axis=1)]

    stable_df.to_csv(outCsv, index=True)
    return stable_df, gene_hits

def repeat_deseq_pseudobulk_analysis_parallel_v2(
    adata,
    ctyps,
    sample_colname,
    cell_type_colname,
    cell_of_interest,
    metacols,
    outCsv,
    N=10,
    min_cells=40,
    max_cells=50,
    min_repeats=5,
    gtf=None,
    n_jobs=4
):
    gene_hits = defaultdict(int)
    final_results = []

    for i in range(N):
        print(f"=== Iteration {i+1}/{N} ===")

        ps_adata = make_pseudobulk_parallel(
            adata,
            sample_colname=sample_colname,
            cell_type_colname=cell_type_colname,
            metacols=metacols,
            min_cells=min_cells,
            max_cells=max_cells,
            random_state=i,
            n_jobs=n_jobs
        )

        # Just in case: check what cell types survived the pseudobulk
        unique_cts = ps_adata.obs[cell_type_colname].unique()
        print(f"Pseudobulk cell types in iteration {i+1}: {unique_cts}")

        # Filter ctyps to those that actually exist after pseudobulk
        valid_ctyps = [ct for ct in ctyps if ct in unique_cts]

        if len(valid_ctyps) == 0:
            print(f"⚠️  No valid cell types to compare in iteration {i+1}. Skipping.")
            continue

        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(run_deseq_for_celltype_v2)(
                ps_adata,
                c,
                cell_type_colname,
                cell_of_interest,
                metacols,
                gtf
            )
            for c in valid_ctyps
        )

        for de in results:
            if de is not None and not de.empty:
                for g, c in zip(de['gene_name'], de['SubType']):
                    gene_hits[(g, c)] += 1
                final_results.append(de)

    if len(final_results) == 0:
        print("⚠️  No DE results found across all iterations. Nothing to save.")
        return None, None

    all_deg_df = pd.concat(final_results)

    stable_genes = [
        (gene, subtype) for (gene, subtype), count in gene_hits.items()
        if count >= min_repeats
    ]

    stable_df = all_deg_df[
        all_deg_df.apply(
            lambda row: (row['gene_name'], row['SubType']) in stable_genes,
            axis=1
        )
    ]

    stable_df.to_csv(outCsv, index=True)
    print(f"✅ Saved stable DEGs to {outCsv}")
    return stable_df, gene_hits

adata = sc.read_h5ad('data/LucaExtended_downloaded.06-21-2024.h5ad')
adata = adata[adata.obs['disease'].isin(['normal', 'lung adenocarcinoma'])]
adata = adata[adata.obs['assay'] == '10x 3\' v2']
gtf_data = parse_gtf('data/Homo_sapiens.GRCh38.104.gtf')

cell_of_interest = "malignant cell"
ctyps = adata.obs['cell_type'].unique()
ctyps = [ct for ct in ctyps if ct != cell_of_interest]
stable_degs, gene_hits_dict = repeat_deseq_pseudobulk_analysis_parallel_v2(
    adata=adata,
    ctyps=ctyps,
    sample_colname='sample',
    cell_type_colname='cell_type',
    cell_of_interest=cell_of_interest,
    metacols=['assay', 'donor_id', 'disease', 'tissue', 'study', 'sex', 'age', 'uicc_stage', 'tumor_stage'],
    N=100,  # number of iterations
    min_cells=40,
    max_cells=50,
    min_repeats=90,
    gtf=gtf_data,  # or None if not available
    n_jobs=10,
    outCsv='results/Luad_PB-DEGs_repeated-N100-minrep90.csv'
)

