import numpy as np
import pandas as pd
import itertools
import scanpy as sc

# # Simulate mock data for 5 genes (in real use, replace with actual data)
# data = pd.DataFrame({
#     "Gene": ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"],
#     "T_mean": [50, 80, 30, 70, 90],            # Tumor expression
#     "Eo": [5, 10, 8, 20, 4],                   # Other lung cell expression
#     "Ee": [3, 6, 4, 15, 2],                    # Epithelial cell expression
#     "B": [0.5, 1.2, 3.0, 0.8, 0.3],            # Brain/CNS expression
#     "C": [4, 3, 2, 5, 6],                      # Differentially expressed cell types
#     "L": [1, 0, 1, 1, 1],                      # Surface localization (1 or 0)
#     "CV_intra": [0.3, 0.5, 0.6, 0.2, 0.4],     # Intra-patient heterogeneity
#     "CV_inter": [0.2, 0.6, 0.4, 0.3, 0.5],     # Inter-patient heterogeneity
#     "U_intra": [0.9, 0.7, 0.6, 0.95, 0.85],    # Uniformity within tumor
#     "U_inter": [0.8, 0.5, 0.4, 0.85, 0.9]      # Uniformity across patients
# })

# Define weight ranges for grid search
# weight_grid = {
#     "alpha": [1],
#     "beta": [1],
#     "gamma": [1],
#     "delta": [2],
#     "lambda_": [2],
#     "eta": [1],
#     "theta": [1],
#     "mu": [2],
#     "nu": [2]
# }

# Create all combinations of weights
# weight_combinations = list(itertools.product(*weight_grid.values()))

# Evaluate each combination and store the results
# results = []
# for weights in weight_combinations:
#     scores = data.apply(lambda row: calculate_cps(row, weights), axis=1)
#     ranked_genes = data.copy()
#     ranked_genes["CPS"] = scores
#     ranked_genes["weights"] = str(weights)
#     ranked_genes_sorted = ranked_genes.sort_values(by="CPS", ascending=False)
#     results.append(ranked_genes_sorted)

# import ace_tools as tools; tools.display_dataframe_to_user(name="CPS Candidate Ranking", dataframe=results[0])



# #sc Scanpy
# def calculate_cps(T_mean, Eo, Ee, C, L, B, CV_intra, CV_inter, U_intra, U_inter,
#                   alpha=1, beta=1, gamma=1, delta=2, lambda_=2, eta=1, theta=1, mu=2, nu=2):
#     return (
#         alpha * np.log2(T_mean / (Eo + 1)) +
#         beta * np.log2(T_mean / (Ee + 1)) +
#         gamma * C +
#         delta * L -
#         lambda_ * np.log2(B + 1) -
#         eta * CV_intra -
#         theta * CV_inter +
#         mu * U_intra +
#         nu * U_inter
#     )



# Load AnnData object
adata = sc.read_h5ad('data/LucaExtended_downloaded.06-21-2024.h5ad')
adata = adata[adata.obs['disease'].isin(['normal', 'lung adenocarcinoma'])]
adata = adata[adata.obs['assay'] == '10x 3\' v2']



# Assume the following exist:
# adata.obs['cell_type'] — cell type labels
# adata.obs['sample'] — patient/sample ID
# adata.raw.X — normalized expression data
# adata.var['gene_name'] — gene symbols

# List of target genes to evaluate
#target_genes = ["MUC1", "CLDN4", "HER2", "GPC3", "EPCAM", "CEACAM5"]  # etc.
# target_genes = [
#     "ENSG00000185499",  # MUC1
#     "ENSG00000189143",  # CLDN4
#     "ENSG00000119888",  # EPCAM
#     "ENSG00000105388",  # CEACAM5
#     "ENSG00000169894",  # MUC3A
#     "ENSG00000131650",  # KREMEN2
#     "ENSG00000165215",  # CLDN3
#     "ENSG00000146648",  # EGFR
#     "ENSG00000062038",  # CDH3
#     "ENSG00000134873",  # CLDN10
#     "ENSG00000130821",  # SLC6A8
#     "ENSG00000204544",  # MUC21
#     "ENSG00000205213",  # LGR4
#     "ENSG00000137975"  # CLCA2
# ]

# Data as list of tuples for correct column naming
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

# Now you can access as:
print(df[["Ensembl_ID", "Gene_Symbol"]])

# Pseudobulk per patient
#grouped = adata.to_df(layer="raw").groupby([adata.obs['sample'], adata.obs['cell_type']]).mean()

# Get raw expression matrix
raw_df = pd.DataFrame(
    adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X,
    index=adata.obs_names,
    columns=adata.raw.var_names)
# Combine with metadata
raw_df["sample"] = adata.obs["sample"].values
raw_df["cell_type"] = adata.obs["cell_type"].values

# Group by and calculate mean
#grouped = raw_df.groupby(["sample", "cell_type"]).mean()


# Aggregate tumor cells only for each patient
tumor_mask = adata.obs['cell_type'].str.contains("malignant cell", case=False)
adata_tumor = adata[tumor_mask, :]


#Get the B
hbca_data = sc.read_h5ad('/mnt/lun3/processed/EbruKocakaya/int_fresh_start/descartes_ebru/Human_Brain_Cell_Atlas_v1.0_all_downloaded.06-22-2024.h5ad')

raw_dfb = pd.DataFrame(
    hbca_data.raw.X.toarray() if hasattr(hbca_data.raw.X, "toarray") else hbca_data.raw.X,
    index=hbca_data.obs_names,
    columns=hbca_data.raw.var_names)
raw_dfb["sample"] = hbca_data.obs["sample"].values
raw_dfb["cell_type"] = hbca_data.obs["cell_type"].values



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


    print(gene, T_mean, CV_intra, CV_inter, U_intra, U_inter)
    gene_stats[gene] = {
        "T_mean": T_mean,
        "CV_intra": CV_intra,
        "CV_inter": CV_inter,
        "U_intra": U_intra,
        "U_inter": U_inter
    }

# Load external data (e.g. from GTEx or Human Protein Atlas)
# external_df = pd.read_csv("external_expression_data.csv")
# external_df should have gene, Eo, Ee, B, L, C columns

degdata = pd.read_csv("results/Luad_PB-DEGs_repeated-allTested.csv")

# Define epithelial cell types
epithelial_labels = {
    "club cell", "multi-ciliated epithelial cell", "type I pneumocyte", "type II pneumocyte"
}

# Label cell types as Ee (epithelial) or Eo (other)
degdata["cell_type_group"] = degdata["SubType"].apply(
    lambda x: "Ee" if x in epithelial_labels else "Eo"
)

# Filter to only your target genes
target_genes = df["Ensembl_ID"].tolist()
deg_filtered = degdata[degdata["Unnamed: 0"].isin(target_genes)]

# Total number of unique cell types tested (you can set this dynamically)
total_celltypes = 30

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

# Display or export
print(result)
# result.to_csv("final_gene_stats_normalized.csv", index=False)
# Optional: result.to_csv("final_gene_stats_with_counts.csv", index=False)

# # Compute average log2FoldChange for each gene in each group
# avg_logfc = (
#     deg_filtered
#     .groupby(["Unnamed: 0", "cell_type_group"])["log2FoldChange"]
#     .mean()
#     .unstack(fill_value=0)  # Fill with 0 if one group is missing
#     .rename_axis("Ensembl_ID")
#     .reset_index()
# )

# # Merge gene symbols
# result = df[["Ensembl_ID", "Gene_Symbol"]].merge(avg_logfc, on="Ensembl_ID", how="left")

# # Optional: rename columns
# result = result.rename(columns={"Ee": "AvgLogFC_Ee", "Eo": "AvgLogFC_Eo"})

# # Display or save
# print(result)
# result.to_csv("avg_logfc_by_group.csv", index=False)


gene_stats_df = pd.DataFrame.from_dict(gene_stats, orient='index').reset_index().rename(columns={"index": "Ensembl_ID"})
final_result = result.merge(gene_stats_df, on="Ensembl_ID", how="left")

final_result["B"] = np.round(np.random.normal(loc=final_result["AvgLogFC_Eo"].mean(), scale=1.0, size=len(final_result)), 6)
final_result["L"] = np.random.choice([0, 1], size=len(final_result))


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

weights = [1, 1, 1, 2, 2, 1, 1, 2, 2]
final_result = calculate_gene_score(final_result, weights)
print(final_result[["Gene_Symbol", "Score"]].sort_values(by="Score", ascending=False))

final_result.to_csv("results/Selected_genes_CPS.csv", index=False)

# external_df = pd.DataFrame({
#     "gene": target_genes,
#     "Eo": [5, 10, 8, 2, 15, 7],
#     "Ee": [3, 6, 4, 1, 9, 6],
#     "B": [0.5, 3.0, 2.5, 0.3, 1.2, 1.0],
#     "L": [1, 1, 0, 1, 0, 1],
#     "C": [5, 4, 3, 6, 2, 3]
# })



# weights = [1, 1, 1, 2, 2, 1, 1, 2, 2]
#weights = weight_combinations
# records = []

# for gene in target_genes:
#     row = {
#         "AvgLogFC_Ee": final_result[final_result["Ensembl_ID"].isin([gene])]["AvgLogFC_Ee"],
#         "AvgLogFC_Eo": final_result[final_result["Ensembl_ID"].isin([gene])]["AvgLogFC_Eo"],
#         "n_upregulated_norm": final_result[final_result["Ensembl_ID"].isin([gene])]["n_upregulated_norm"],
#         "L": final_result[final_result["Ensembl_ID"].isin([gene])]["L"],
#         "B": final_result[final_result["Ensembl_ID"].isin([gene])]["B"],
#         "CV_intra": final_result[final_result["Ensembl_ID"].isin([gene])]["CV_intra"],
#         "CV_inter": final_result[final_result["Ensembl_ID"].isin([gene])]["CV_inter"],
#         "U_intra": final_result[final_result["Ensembl_ID"].isin([gene])]["U_intra"],
#         "U_inter": final_result[final_result["Ensembl_ID"].isin([gene])]["U_inter"]
#     }
#     print(row)
#     cps = calculate_cps(row, weights)
#     records.append({
#         "Gene": gene,
#         "CPS": cps,
#         **stats,
#         **ext.to_dict()
#     })

# Merge all features
# records = []
# for gene in target_genes:
#     #stats = gene_stats[gene]
#     #ext = external_df[external_df["gene"] == gene].iloc[0]
#     cps = calculate_cps(
#         AvgLogFC_Ee=final_result["AvgLogFC_Ee"],
#         AvgLogFC_Eo=final_result["AvgLogFC_Eo"],
#         C=final_result["n_upregulated_norm"], # This is C
#         L=final_result["L"],
#         B=final_result["B"],
#         CV_intra=final_result["CV_intra"],
#         CV_inter=final_result["CV_inter"],
#         U_intra=final_result["U_intra"],
#         U_inter=final_result["U_inter"]
#     )
#     records.append({
#         "Gene": gene,
#         "CPS": cps,
#         **stats,
#         **ext.to_dict()
#     })

#from scripts.subroutines import calculate_gene_score

#result_df = pd.DataFrame(records).sort_values(by="CPS", ascending=False)
