"""
===========================================================================
Pseudobulk Generator
===========================================================================

This script generates a pseudobulk gene expression matrix from single-cell data stored in an AnnData (.h5ad) object.

Functionality:
1. Groups cells by a specified patient/sample column in `adata.obs`.
2. Sums gene expression counts across all cells for each patient/sample.
3. Outputs a CSV file where rows are patients/samples and columns are genes.

Notes:
- Supports sparse matrices (csr_matrix) efficiently to handle large datasets.
- If the expression matrix is dense, it will warn the user (may be memory intensive).
- Be aware for Statescope you need to transpose the output!

Inputs (via command-line arguments):
--input_h5ad     Path to input .h5ad file
--patient_col    Column name in adata.obs used to group cells (e.g., 'case_id')
--output_csv     Path to save the resulting pseudobulk CSV

Outputs:
- CSV file containing summed gene expression per patient/sample

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)

"""

import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
import argparse


def create_pseudobulk(adata, patient_id, pseudobulk_file):
    print('Creating pseudobulk...')
    # Group by 'case_id' and sum the gene expression counts
    print('Saving adata.X into expression_matrix variable...')
    expression_matrix = adata.X

    # Check if matrix is sparse (if not and the file is too large the script might crash)
    print('Checking if sparse matrix...')
    if isinstance(expression_matrix, csr_matrix):
        # Create a DataFrame where rows are 'case_id' and columns are gene names, sum each 
        print('Creating data frame...')
        summed_expression = pd.DataFrame.sparse.from_spmatrix(expression_matrix, index=adata.obs[patient_id], columns=adata.var_names).groupby(level=0).sum()
        summed_expression.index.name = None
    else:
        # If it's not sparse (should not happen if adata.X is large), use regular dense handling
        print("Expression matrix is not sparse.")

    # Print resulting DataFrame
    print(summed_expression.head())
    print(summed_expression.index)

    # Save
    print('Saving to csv...')
    summed_expression.to_csv(pseudobulk_file)
    return summed_expression

# -------------------------
# Command-line interface
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input_h5ad", required=True)
parser.add_argument("--patient_col", required=True)
parser.add_argument("--output_csv", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input_h5ad)
create_pseudobulk(adata, args.patient_col, args.output_csv)
