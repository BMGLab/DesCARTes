import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

# warnings.filterwarnings("ignore")

def add_brain_expression_column(final_df, brain_expr_csv, b_column_name="B"):
    """
    Add brain expression values as column B for CPS scoring.

    Parameters:
    ----------
    final_df : pd.DataFrame
        The main DataFrame containing Ensembl_IDs and gene information.
    brain_expr_csv : str
        Path to the CSV file containing brain average expression with columns ['Gene', 'AverageExpression'].
    b_column_name : str
        Name of the output column in final_df to store the brain expression (default: 'B').

    Returns:
    -------
    pd.DataFrame
        Updated DataFrame with a new column (default 'B') containing brain expression values.
    """
    brain_expr = pd.read_csv(brain_expr_csv)
    brain_expr = brain_expr.rename(columns={
        "Gene": "Ensembl_ID",
        "AverageExpression": "BrainExpression"
    })
    merged_df = final_df.merge(brain_expr, on="Ensembl_ID", how="left")
    merged_df[b_column_name] = merged_df["BrainExpression"].fillna(0)
    return merged_df

def add_surface_annotation_column(final_df, surfaceome_path, l_column_name="L"):
    """
    Annotate genes with a binary surface marker label based on a surfaceome reference list.

    Parameters:
    ----------
    final_df : pd.DataFrame
        The main DataFrame containing an 'Ensembl_ID' column.
    surfaceome_path : str
        Path to the surfaceome annotation file (e.g., table_S3_surfaceome.txt).
    l_column_name : str
        Name of the output column in final_df to store the surfaceome label (default: 'L').

    Returns:
    -------
    pd.DataFrame
        Updated DataFrame with a new binary column indicating surfaceome membership.
    """
    surface_df = pd.read_csv(surfaceome_path, sep="\t", low_memory=False)
    surface_ensembl_ids = set(surface_df["Ensembl gene"].dropna().unique())
    final_df[l_column_name] = final_df["Ensembl_ID"].apply(lambda eid: 1 if eid in surface_ensembl_ids else 0)
    return final_df

# === Example usage ===
# filter_luad_and_normal("/path/to/full_file.h5ad")
def filter_luad_and_normal(
    input_h5ad: str,
    output_h5ad: str = "data/adata_LUAD_normal_filtered.h5ad",
    luad_keywords = ["lung adenocarcinoma", "LUAD"],
    normal_keywords = ["normal", "healthy"]
):
    """
    Filters LUAD and normal samples from an .h5ad file and saves a new file.

    Parameters:
    - input_h5ad (str): Path to the full input AnnData file.
    - output_h5ad (str): Path where the filtered AnnData will be saved.
    - luad_keywords (list): List of strings representing LUAD labels.
    - normal_keywords (list): List of strings representing normal/healthy labels.
    """
    import anndata as ad
    import os
    if os.path.exists(output_h5ad):
        print(f"\tFiltered file already exists at: {output_h5ad}")
        print("\tSkipping filtering step.")
        return

    print("\tLoading full AnnData...")
    adata = sc.read_h5ad(input_h5ad)

    # Identify the appropriate .obs column
    obs_col = "condition" if "condition" in adata.obs.columns else "disease"
    print(f"\tFiltering cells using obs column: '{obs_col}'")

    # Create a filter mask using case-insensitive match
    keywords = luad_keywords + normal_keywords
    mask = adata.obs[obs_col].str.lower().str.contains('|'.join(map(str.lower, keywords)), na=False)
    adata_filtered = adata[mask].copy()

    print(f"\tCells before filtering: {adata.n_obs}")
    print(f"\tCells after filtering:  {adata_filtered.n_obs}")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)

    print(f"\tSaving filtered AnnData to: {output_h5ad}")
    adata_filtered.write(output_h5ad)
    print("\tDone. Filtered file ready for downstream analysis.")


def calculate_ubiquity_scores(
    adata,
    genes,
    threshold=1.0,
    tumor_identifier="malignant cell",
    sample_col="sample",
    cell_type_col="cell_type",
    expr_df=None,
    q2_masks=None
):
    """
    Calculate U_intra and U_inter for a list of genes, optionally using Q2-filtered cells.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (must include raw)
    genes : list
        List of gene identifiers (Ensembl IDs)
    threshold : float, optional
        Expression threshold (default=1.0 for raw counts, or 0.0 for z-scores)
    tumor_identifier : str, optional
        String to identify tumor cells (default="malignant cell")
    sample_col : str, optional
        Column in .obs indicating patient/sample ID
    cell_type_col : str, optional
        Column in .obs with cell type annotations
    expr_df : pd.DataFrame, optional
        Precomputed expression DataFrame (raw or z-scored), indexed by cell, columns = genes
    q2_masks : dict, optional
        Optional dictionary: gene → Boolean mask (True = include cell in Q2 calc)

    Returns
    -------
    (dict, dict)
        Tuple of (U_intra_scores, U_inter_scores)
    """

    # Filter tumor cells
    tumor_mask = adata.obs[cell_type_col].str.contains(tumor_identifier, case=False)
    adata_tumor = adata[tumor_mask].copy()
    if expr_df is None:
        expr_df = pd.DataFrame(
            adata_tumor.raw.X.toarray() if hasattr(adata_tumor.raw.X, "toarray")
            else adata_tumor.raw.X,
            index=adata_tumor.obs_names,
            columns=adata_tumor.raw.var_names
        )
    # Add metadata
    expr_df["sample"] = adata_tumor.obs[sample_col].values
    # Output containers
    u_intra_scores = {}
    u_inter_scores = {}
    for gene in genes:
        try:
            gene_expr = expr_df[gene]
        except KeyError:
            # Skip missing gene
            u_intra_scores[gene] = np.nan
            u_inter_scores[gene] = np.nan
            continue
        # Optional Q2 filter
        if q2_masks and gene in q2_masks:
            gene_mask = q2_masks[gene]
            gene_expr = gene_expr[gene_mask]
            sample_ids = expr_df.loc[gene_mask, "sample"]
        else:
            sample_ids = expr_df["sample"]

        # Group gene expression by sample
        expr_by_patient = gene_expr.groupby(sample_ids, observed=False)

        # Per-patient coverage: fraction of cells expressing the gene (expression > threshold)
        patient_coverage = expr_by_patient.apply(lambda x: (x > threshold).sum() / len())

        # Percentage of expressing cells per patient, then average
        u_intra_scores[gene] = 100 * patient_coverage.mean()

        # Percentage of all tumor cells expressing the gene
        u_inter_scores[gene] = 100 * (gene_expr > threshold).sum() /len()
    return u_intra_scores, u_inter_scores

def calculate_gene_score(df, weights):
    alpha, beta, gamma, delta, lambda_, eta, theta, mu, nu = weights
    B_safe = df["B"].clip(lower=1e-6)
    df["Score"] = (delta * df["L"]) * (
        alpha * np.log2(df["T_mean"] / (df["AvgLogFC_Eo"] + 1)) +
        beta * np.log2(df["T_mean"] / (df["AvgLogFC_Ee"] + 1)) +
        gamma * df["n_upregulated_norm"] -
        lambda_ * np.log2(B_safe + 1) -
        eta * df["CV_intra_Q2"] -
        theta * df["CV_inter_Q2"] +
        mu * df["U_intra_Q2"] +
        nu * df["U_inter_Q2"]
    )
    return df

# Explore correlation between variability and ubiquity for genes — for example:
#   - Genes with high ubiquity & low variability could be “housekeeping genes.”
#   - Genes with low ubiquity & high variability could be more specialized.
def plot_cv_vs_ubiquity(df, x_col="CV_intra", y_col="U_intra", outdir="results/plots"):
    os.makedirs(outdir, exist_ok=True)
    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=df, x=x_col, y=y_col, alpha=0.7, edgecolor=None)
    plt.title(f"{x_col} vs. {y_col}")
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.tight_layout()
    fname = f"{x_col}_vs_{y_col}.png"
    plt.savefig(os.path.join(outdir, fname), dpi=150)
    plt.close()

# - See how filtering changes the spread and shape of a metric.
# - Check if filtering skews values toward higher/lower CV or ubiquity.
def plot_raw_vs_q2_distributions(df, metric_base, title_suffix="", outdir="results/plots"):
    os.makedirs(outdir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    sns.histplot(df[f"{metric_base}"], bins=50, kde=True, ax=axes[0], color='steelblue')
    axes[0].set_title(f"Raw {metric_base} {title_suffix}")
    axes[0].set_xlabel(metric_base)
    sns.histplot(df[f"{metric_base}_Q2"], bins=50, kde=True, ax=axes[1], color='darkorange')
    axes[1].set_title(f"Q2-filtered {metric_base} {title_suffix}")
    axes[1].set_xlabel(f"{metric_base}_Q2")
    plt.tight_layout()
    fname = f"{metric_base}_raw_vs_Q2{title_suffix.replace(' ', '_')}.png"
    plt.savefig(os.path.join(outdir, fname), dpi=150)
    plt.close(fig)

# Directly compare shifts in distribution between raw and filtered data without having to look at two separate subplots.
def overlay_density(df, metric_base, label1="Raw", label2="Q2-filtered", outdir="results/plots"):
    os.makedirs(outdir, exist_ok=True)
    plt.figure(figsize=(6, 4))
    sns.kdeplot(df[f"{metric_base}"], label=label1, linewidth=2)
    sns.kdeplot(df[f"{metric_base}_Q2"], label=label2, linewidth=2, linestyle='--')
    plt.title(f"{metric_base} vs. {metric_base}_Q2")
    plt.xlabel(metric_base)
    plt.legend()
    plt.tight_layout()
    fname = f"{metric_base}_overlay_density.png"
    plt.savefig(os.path.join(outdir, fname), dpi=150)
    plt.close()

def main():
    # --- Parameters ---
    weights = [2, 2, 3, 1, 1, 1, 1, 1, 1]
    epithelial_labels = {
        "club cell", "multi-ciliated epithelial cell", "type I pneumocyte", "type II pneumocyte"
    }
    total_celltypes = 30

    # --- Data Loading ---
    print("Loading AnnData...")
    filter_luad_and_normal('data/LucaExtended_downloaded.06-21-2024.h5ad')
    adata = sc.read_h5ad('data/adata_LUAD_normal_filtered.h5ad')
    adata = adata[adata.obs['disease'].isin(['normal', 'lung adenocarcinoma'])]
    adata = adata[adata.obs['assay'] == '10x 3\' v2']
    print(f"\tNumber of cells: {adata.n_obs}, Number of genes: {adata.n_vars}")

    print("Loading pseudo-bulk results...")
    degdata = pd.read_csv("results/Luad_PB-DEGs_repeated-allTested.csv")
    print(f"\tNumber of DEGs: {degdata.shape[0]}")

    # --- Gene Symbol Mapping ---
    all_genes = degdata["Unnamed: 0"].unique()
    df = pd.DataFrame({"Ensembl_ID": all_genes})
    if "gene_name" in degdata.columns:
        id_to_symbol = dict(zip(degdata["Unnamed: 0"], degdata["gene_name"]))
        df["Gene_Symbol"] = df["Ensembl_ID"].map(id_to_symbol)
    elif "gene_symbols" in adata.raw.var.columns:
        id_to_symbol = dict(zip(adata.raw.var_names, adata.raw.var["gene_symbols"]))
        df["Gene_Symbol"] = df["Ensembl_ID"].map(id_to_symbol)
    else:
        df["Gene_Symbol"] = df["Ensembl_ID"]

    target_genes = df["Ensembl_ID"].tolist()

    # --- Extract Tumor Cells ---
    tumor_mask = adata.obs['cell_type'].str.contains("malignant cell", case=False)
    adata_tumor = adata[tumor_mask, :]
    print(f"\tNumber of tumor cells: {adata_tumor.n_obs}, Number of genes: {adata_tumor.n_vars}")

    tumor_expr_df = adata_tumor.to_df()

    # --- Z-Score Normalization ---
    print("Applying z-score normalization...")
    tumor_expr_z = tumor_expr_df.apply(zscore, axis=0, nan_policy='omit')

    print("Generating Q1–Q3 masks for each gene...")
    q2_masks = {
        gene: (tumor_expr_z[gene] >= tumor_expr_z[gene].quantile(0.25)) &
              (tumor_expr_z[gene] <= tumor_expr_z[gene].quantile(0.75))
        for gene in tumor_expr_z.columns
    }

    # --- Modified Ubiquity Scores ---
    print("Calculating Q2-filtered ubiquity scores...")
    u_intra_q2, u_inter_q2 = calculate_ubiquity_scores(
        adata=adata_tumor,
        genes=target_genes,
        threshold=0,  # z-score space
        u_inter_threshold=0.5,
        tumor_identifier="malignant cell",
        sample_col="sample",
        cell_type_col="cell_type"
    )

    print("Calculating raw ubiquity scores...")
    u_intra_raw, u_inter_raw = calculate_ubiquity_scores(
        adata=adata_tumor,
        genes=target_genes,
        threshold=1.0  # raw counts
    )

    # --- Gene-level Stats: T, CV, U ---
    print("Calculating per-gene stats...")
    gene_stats = {}
    for gene in target_genes:
        full_expr = tumor_expr_df[gene]         # Raw expression values across all tumor cells
        z_expr = tumor_expr_z[gene]             # Z-score normalized expression for the gene
        mask = q2_masks[gene]                   # Boolean mask to keep cells in the IQR (Q2 filter)

        # Group raw expression by sample (patient-level)
        full_grouped = full_expr.groupby(adata_tumor.obs["sample"], observed=False)
        # Group Q2-filtered z-scored expression by sample
        q2_grouped = z_expr[mask].groupby(adata_tumor.obs.loc[mask, "sample"], observed=False)

        # === Raw Expression-Based Statistics ===
        # T_mean: Average expression of gene across patients (mean of per-sample means)
        T_mean = full_grouped.apply(np.mean).mean()
        # CV_intra: Mean of within-patient CVs (per sample)
        CV_intra = full_grouped.apply(lambda x: np.std(x) / (np.mean(x) + 1e-6)).mean()
        # CV_inter: CV of mean expression values across patients (between-patient variability)
        CV_inter = np.std(full_grouped.apply(np.mean)) / (np.mean(full_grouped.apply(np.mean)) + 1e-6)

        # === Q2-Filtered Z-Score-Based Statistics ===
        # CV_intra_Q2: Mean of CVs across patients, but calculated on Q2-filtered z-scores
        CV_intra_Q2 = (
            q2_grouped.apply(lambda x: np.std(x) / (np.mean(x) + 1e-6)).mean()
            if len(q2_grouped) > 0 else np.nan
        )
        # CV_inter_Q2: CV of the sample-level means (after Q2 filtering), captures between-patient consistency
        CV_inter_Q2 = (
            np.std(q2_grouped.apply(np.mean)) /
            (np.mean(q2_grouped.apply(np.mean)) + 1e-6)
            if len(q2_grouped) > 0 else np.nan
        )
        # Store all calculated metrics for the gene
        gene_stats[gene] = {
            "T_mean": T_mean,
            "CV_intra": CV_intra,
            "CV_inter": CV_inter,
            "CV_intra_Q2": CV_intra_Q2,
            "CV_inter_Q2": CV_inter_Q2,
            "U_intra": u_intra_raw[gene],
            "U_inter": u_inter_raw[gene],
            "U_intra_Q2": u_intra_q2[gene],
            "U_inter_Q2": u_inter_q2[gene],
        }

    # --- DESeq2 Processing ---
    degdata["cell_type_group"] = degdata["SubType"].apply(
        lambda x: "Ee" if x in epithelial_labels else "Eo")
    deg_filtered = degdata[degdata["Unnamed: 0"].isin(target_genes)]

    print("Computing Avg LogFC...")
    avg_logfc = (
        deg_filtered
        .groupby(["Unnamed: 0", "cell_type_group"])["log2FoldChange"]
        .mean()
        .unstack(fill_value=0)
        .rename_axis("Ensembl_ID")
        .reset_index()
    )

    print("Counting # of upregulated SubTypes per gene")
    sig_up = deg_filtered[
        (deg_filtered["padj"] < 0.05) &
        (deg_filtered["log2FoldChange"] > 0)
    ]
    upreg_count = (
        sig_up.groupby("Unnamed: 0")["SubType"]
        .nunique()
        .reset_index()
        .rename(columns={"Unnamed: 0": "Ensembl_ID", "SubType": "n_upregulated_celltypes"})
    )
    upreg_count["n_upregulated_norm"] = upreg_count["n_upregulated_celltypes"] / total_celltypes

    # --- Merge & Compute Score ---
    print("Merging data and computing CPS...")
    result = df[["Ensembl_ID", "Gene_Symbol"]] \
        .merge(avg_logfc, on="Ensembl_ID", how="left") \
        .rename(columns={"Ee": "AvgLogFC_Ee", "Eo": "AvgLogFC_Eo"}) \
        .merge(upreg_count, on="Ensembl_ID", how="left")

    result["n_upregulated_celltypes"] = result["n_upregulated_celltypes"].fillna(0).astype(int)
    result["n_upregulated_norm"] = result["n_upregulated_norm"].fillna(0.0)

    gene_stats_df = pd.DataFrame.from_dict(gene_stats, orient='index').reset_index().rename(columns={"index": "Ensembl_ID"})

    plot_cv_vs_ubiquity(gene_stats_df, "CV_inter_Q2", "U_inter_Q2")
    plot_cv_vs_ubiquity(gene_stats_df, "CV_intra_Q2", "U_intra_Q2")  # Q2 version

    # --- Plot distributions for CVs and Us ---
    print("Plotting distributions for CVs and Us (raw vs Q2-filtered)...")
    plot_raw_vs_q2_distributions(gene_stats_df, "CV_intra", "(within-patient)")
    plot_raw_vs_q2_distributions(gene_stats_df, "CV_inter", "(between-patient)")
    plot_raw_vs_q2_distributions(gene_stats_df, "U_intra", "(within-patient ubiquity)")
    plot_raw_vs_q2_distributions(gene_stats_df, "U_inter", "(between-patient ubiquity)")

    # Optional: Overlay density plots for direct comparison
    overlay_density(gene_stats_df, "CV_intra")
    overlay_density(gene_stats_df, "CV_inter")
    overlay_density(gene_stats_df, "U_intra")
    overlay_density(gene_stats_df, "U_inter")

    final_result = result.merge(gene_stats_df, on="Ensembl_ID", how="left")
    final_result = add_brain_expression_column(final_result, "data/average_Brain.gene_expression.csv")
    final_result = add_surface_annotation_column(final_result, "data/table_S3_surfaceome.txt")
    final_result = calculate_gene_score(final_result, weights)

    final_result.to_csv("results/Selected_genes_CPS_with_Q2.csv", index=False)
    print("✅ CPS scores calculated and saved.")

if __name__ == "__main__":
    main()