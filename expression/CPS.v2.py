import numpy as np
import pandas as pd
import scanpy as sc
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

def calculate_gene_score(df, weights):
    alpha, beta, gamma, delta, lambda_, eta, theta, mu, nu = weights
    B_safe = df["B"].clip(lower=1e-6)
    df["Score"] = (
        alpha * df["AvgLogFC_Eo"] +
        beta * df["AvgLogFC_Ee"] +
        gamma * df["n_upregulated_norm"] +
        delta * df["L"] -
        lambda_ * np.log2(B_safe + 1) -
        eta * df["CV_intra"] -
        theta * df["CV_inter"] +
        mu * df["U_intra"] +
        nu * df["U_inter"]
    )
    return df

# define the weights for the CPS calculation. !!! randomly generated for demonstration purposes
#[alpha, beta, gamma, delta, lambda_, eta, theta, mu, nu]
weights = [2, 2, 3, 1, 3, 2, 2, 1, 1]

# Load AnnData object
print("Loading AnnData...")
adata = sc.read_h5ad('data/LucaExtended_downloaded.06-21-2024.h5ad')
adata = adata[adata.obs['disease'].isin(['normal', 'lung adenocarcinoma'])]
adata = adata[adata.obs['assay'] == '10x 3\' v2']
print(f"Number of cells: {adata.n_obs}, Number of genes: {adata.n_vars}")

# Load external data (pseudo-bulk DESeq2 results)
print("Loading pseudo-bulk results...")
degdata = pd.read_csv("results/Luad_PB-DEGs_repeated-allTested.csv")
print(f"Number of DEGs: {degdata.shape[0]}\n")

# Define epithelial cell types
epithelial_labels = {
    "club cell", "multi-ciliated epithelial cell", "type I pneumocyte", "type II pneumocyte"
}

# Total number of unique cell types tested in DESeq2 analysis
total_celltypes = 30 # Adjust this based on your dataset or findout dynamically

# Data as list of tuples for correct column naming --> change to all genes of interest
data = [
    ("ENSG00000185499", "MUC1"),
    ("ENSG00000189143", "CLDN4"),
    ("ENSG00000119888", "EPCAM"),
    ("ENSG00000105388", "CEACAM5"),
    ("ENSG00000169894", "MUC3A"),
    ("ENSG00000131650", "KREMEN2"),
    ("ENSG00000165215", "CLDN3"),
    ("ENSG00000146648", "EGFR"),
    ("ENSG00000062038", "CDH3"),
    ("ENSG00000134873", "CLDN10"),
    ("ENSG00000130821", "SLC6A8"),
    ("ENSG00000204544", "MUC21"),
    ("ENSG00000205213", "LGR4"),
    ("ENSG00000137975", "CLCA2")
]

# Create DataFrame with correct column names
df = pd.DataFrame(data, columns=["Ensembl_ID", "Gene_Symbol"])
target_genes = df["Ensembl_ID"].tolist()

print("-------------------------------------------------------")
print(df)
print("-------------------------------------------------------")

# Get raw expression matrix
raw_df = pd.DataFrame(
    adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X,
    index=adata.obs_names,
    columns=adata.raw.var_names)
# Combine with metadata
raw_df["sample"] = adata.obs["sample"].values
raw_df["cell_type"] = adata.obs["cell_type"].values

# Aggregate tumor cells only for each patient
tumor_mask = adata.obs['cell_type'].str.contains("malignant cell", case=False)
adata_tumor = adata[tumor_mask, :]
print(f"Number of tumor cells: {adata_tumor.n_obs}, Number of genes: {adata_tumor.n_vars}")

# Calculate per-gene tumor expression metrics
gene_stats = {}
for gene in target_genes:
    gene_expr_by_patient = adata_tumor.to_df()[gene].groupby(adata_tumor.obs['sample'])
    expr_matrix = gene_expr_by_patient.apply(lambda x: x.values)
    T_mean = expr_matrix.apply(np.mean).mean()
    CV_intra = expr_matrix.apply(lambda x: np.std(x) / np.mean(x + 1e-6)).mean()
    CV_inter = np.std(expr_matrix.apply(np.mean)) / (np.mean(expr_matrix.apply(np.mean)) + 1e-6)
    U_intra = expr_matrix.apply(lambda x: np.mean(x > 1)).mean()
    U_inter = np.mean(expr_matrix.apply(lambda x: np.mean(x > 1)) > 0.5)
    gene_stats[gene] = {
        "T_mean": T_mean,
        "CV_intra": CV_intra,
        "CV_inter": CV_inter,
        "U_intra": U_intra,
        "U_inter": U_inter
    }

# Label cell types as Ee (epithelial) or Eo (other)
degdata["cell_type_group"] = degdata["SubType"].apply(
    lambda x: "Ee" if x in epithelial_labels else "Eo")

# Filter to only your target genes
deg_filtered = degdata[degdata["Unnamed: 0"].isin(target_genes)]

print(f"Number of DEGs for which CPS will be calculated: {deg_filtered.shape[0]}")


# --- Part 1: Compute Avg LogFC ---
avg_logfc = (
    deg_filtered
    .groupby(["Unnamed: 0", "cell_type_group"])["log2FoldChange"]
    .mean()
    .unstack(fill_value=0)
    .rename_axis("Ensembl_ID")
    .reset_index()
)

# --- Part 2: Count # of upregulated SubTypes per gene ---
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

# Add normalized count column
upreg_count["n_upregulated_norm"] = upreg_count["n_upregulated_celltypes"] / total_celltypes

# --- Part 3: Merge everything ---
result = df[["Ensembl_ID", "Gene_Symbol"]] \
    .merge(avg_logfc, on="Ensembl_ID", how="left") \
    .rename(columns={"Ee": "AvgLogFC_Ee", "Eo": "AvgLogFC_Eo"}) \
    .merge(upreg_count, on="Ensembl_ID", how="left")

# Fill missing values
result["n_upregulated_celltypes"] = result["n_upregulated_celltypes"].fillna(0).astype(int)
result["n_upregulated_norm"] = result["n_upregulated_norm"].fillna(0.0)

# print(result)
# Add gene statistics
gene_stats_df = pd.DataFrame.from_dict(gene_stats, orient='index').reset_index().rename(columns={"index": "Ensembl_ID"})
final_result = result.merge(gene_stats_df, on="Ensembl_ID", how="left")

final_result = add_brain_expression_column(final_result, "data/average_Brain.gene_expression.csv")
final_result = add_surface_annotation_column(final_result, "data/table_S3_surfaceome.txt")

final_result = calculate_gene_score(final_result, weights)
print(final_result[["Gene_Symbol", "Score"]].sort_values(by="Score", ascending=False))

final_result.to_csv("results/Selected_genes_CPS_t2.csv", index=False)
print("CPS scores calculated and saved to results/Selected_genes_CPS_t2.csv")
