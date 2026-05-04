"""
Script: TCGA Bulk RNA-seq Sample Filtering

Description:
This script loads a TCGA bulk RNA-seq transcriptome matrix, removes specified samples 
(e.g., normal samples or those without primary tumor data), and overwrites the original 
file with the filtered dataset. It also reports the number of samples before and after filtering.

Author: Mabel Pronk
"""

import pandas as pd

# 1. Load the bulk RNA-seq data matrix (genes x samples)
file_path = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/input/transcriptome_matrix_TCGA_clean_column.csv"
df = pd.read_csv(file_path, index_col=0)

# 2. Define samples to exclude
# These samples correspond to normal tissue or lack primary tumor data
exclude_samples = [
    'TCGA-76-4927', 
    'TCGA-06-5417', 
    'TCGA-06-0167', 
    'TCGA-06-0139', 
    'TCGA-06-0178'
]

# 3. Filter the dataframe
# Retain only samples NOT listed in exclude_samples
original_count = df.shape[1]
df_filtered = df.drop(columns=[s for s in exclude_samples if s in df.columns])

# Print summary of filtering
print(f"Original sample count: {original_count}")
print(f"Final sample count: {df_filtered.shape[1]}")
print(f"Removed: {[s for s in exclude_samples if s in df.columns]}")

# 4. Save the filtered data (overwrite original file)
df_filtered.to_csv(file_path)
print("File updated and saved successfully.")
