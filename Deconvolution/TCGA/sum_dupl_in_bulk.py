"""
Script: Summation of Duplicate Gene Entries in Bulk RNA-seq Data

Description:
This script aggregates duplicated gene entries in a TCGA bulk RNA-seq expression matrix 
by summing their expression values across all samples. Prior to running this script, 
it is recommended to verify whether duplicated genes overlap with signature marker genes 
(e.g., using 'check_dupl_bulk_in_signature.py'), as summing marker genes may affect downstream analyses.

The script groups rows by gene name (index), collapses duplicates by summation, and 
overwrites the original file with a cleaned matrix containing unique gene identifiers.

Author: Mabel pronk
"""

# First use script: check_dupl_bulk_in_signature.py to verify whether duplicates are present 
# in the signature matrix. If no critical conflicts are found, proceed with this script.

import pandas as pd

# 1. Load the bulk RNA-seq expression matrix (genes x samples)
file_path = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/input/transcriptome_matrix_TCGA_clean_column.csv"
df = pd.read_csv(file_path, index_col=0)

print(f"Original shape (with duplicates): {df.shape}")

# 2. Aggregate duplicated gene entries
# Group by gene names (index) and sum expression values across all samples
# This collapses multiple rows of the same gene into a single row
df_summed = df.groupby(df.index).sum()

print(f"New shape (unique genes): {df_summed.shape}")
print(f"Number of duplicate rows collapsed: {len(df) - len(df_summed)}")

# 3. Save the cleaned dataset (overwrite original file)
df_summed.to_csv(file_path)

print("Success: Bulk data is now summed and contains unique gene identifiers. You can now restart your Statescope run.")
