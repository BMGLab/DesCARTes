import scanpy as sc
import pandas as pd
import numpy as np

# Load data (do this only once!)
hbca_data = sc.read_h5ad('/mnt/lun3/processed/EbruKocakaya/int_fresh_start/descartes_ebru/Human_Brain_Cell_Atlas_v1.0_all_downloaded.06-22-2024.h5ad')

# Choose the expression matrix: `.X` or `hbca_data.raw.X` depending on your pipeline
# Assuming normalized values are in .X
X = hbca_data.X
gene_names = hbca_data.var_names

# Compute average expression across all cells
if not isinstance(X, np.ndarray):  # likely sparse
    avg_expr = X.mean(axis=0).A1  # convert to 1D array
else:
    avg_expr = X.mean(axis=0)

# Create DataFrame: rows = genes, column = average expression
avg_expr_df = pd.DataFrame({
    "Gene": gene_names,
    "AverageExpression": avg_expr
})

print(avg_expr_df.head())  # Display first few rows for verification

# Save to disk for future use
avg_expr_df.to_csv("data/average_Brain.gene_expression.csv", index=False)
print("Average gene expression computed and saved to data/average_Brain.gene_expression.csv")


#Do the same for Tabula sapiens
tabula_sapiens_data = sc.read_h5ad('/mnt/lun3/processed/EbruKocakaya/int_fresh_start/descartes_ebru/TabulaSapiens_all_downloaded.06-21-2024.h5ad')

X = tabula_sapiens_data.X
gene_names = tabula_sapiens_data.var_names

# Compute average expression across all cells
if not isinstance(X, np.ndarray):  # likely sparse
    avg_expr = X.mean(axis=0).A1  # convert to 1D array
else:
    avg_expr = X.mean(axis=0)

# Create DataFrame: rows = genes, column = average expression
avg_expr_df = pd.DataFrame({
    "Gene": gene_names,
    "AverageExpression": avg_expr
})

print(avg_expr_df.head())  # Display first few rows for verification

# Save to disk for future use
avg_expr_df.to_csv("data/average_TabulaSapiens.gene_expression.csv", index=False)
print("Average gene expression computed and saved to data/average_TabulaSapiens.gene_expression.csv")

