# ============================================================
# Script: Construct AnnData Objects from Raw 10x Matrix Files
# ============================================================
#
# Description:
# This script processes raw single-cell RNA-seq data stored in 10x Genomics
# Matrix Market format (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
# and converts each sample into a structured AnnData (.h5ad) object.
#
# Each sample is processed individually and saved as a separate file.
# The script also performs basic preprocessing steps such as gene duplicate
# handling and quality control metric computation.
#
# Workflow:
# 1. Identify all sample directories in the input dataset.
# 2. Locate and load the raw 10x Genomics matrix files for each sample.
# 3. Construct an AnnData object from the sparse expression matrix.
# 4. Extract gene names and cell barcodes.
# 5. Detect and export duplicated gene symbols.
# 6. Ensure gene names are unique within each sample.
# 7. Compute basic QC metrics (e.g., mitochondrial gene percentage).
# 8. Add sample-level metadata to AnnData observations.
# 9. Standardize observation names for uniqueness across samples.
# 10. Save each processed sample as an individual .h5ad file.
#
# Inputs:
# - matrix.mtx.gz (expression matrix)
# - barcodes.tsv.gz (cell identifiers)
# - features.tsv.gz (gene identifiers)
#
# Outputs:
# - One AnnData (.h5ad) file per sample
# - Optional CSV files listing duplicated gene names per sample
#
# ============================================================

import os
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from anndata import AnnData

#-------------------------------------------------------------------------------
# 0. Prepare variables
#-------------------------------------------------------------------------------

study_id = "_GSE222522"
cohort = 'GSE222522'
#diagnosis = 'Glioblastoma'
#idh = 'WT'
#pt = 'Primary'
# List all samples

base_data_dir = "/net/beegfs/users/P086608/scRNA_glioma/GSE222522/extracted"
samples = sorted(os.listdir(base_data_dir))

output_dir = "/net/beegfs/users/P086608/scRNA_glioma/GSE222522/adata_files"
os.makedirs(output_dir, exist_ok=True)

#-------------------------------------------------------------------------------
# 1. Loop over samples
#-------------------------------------------------------------------------------

for sample in samples:
    print(f"\nProcessing sample: {sample}")

    # Detect inner folder (NGB1, IMP3, etc.)
    subfolders = [f for f in os.listdir(os.path.join(base_data_dir, sample)) 
                  if os.path.isdir(os.path.join(base_data_dir, sample, f))]
    if len(subfolders) != 1:
        print(f"⚠️ Warning: {sample} has {len(subfolders)} subfolders, skipping.")
        continue
    subfolder = subfolders[0]

    data_dir = os.path.join(base_data_dir, sample, subfolder, "filtered_feature_bc_matrix")
    output_file = os.path.join(output_dir, f"{cohort}_{subfolder}.h5ad")

    matrix_path = os.path.join(data_dir, 'matrix.mtx.gz') 
    barcodes_path = os.path.join(data_dir, 'barcodes.tsv.gz') 
    features_path = os.path.join(data_dir, 'features.tsv.gz') 

    # Load matrix
    print("Loading matrix...")
    matrix = mmread(matrix_path).tocsc().T

    # Load metadata
    print("Loading genes and barcodes...")
    gene_table = pd.read_csv(features_path, sep="\t", header=None)
    barcodes = pd.read_csv(barcodes_path, sep='\t', header=None)[0].values
    gene_symbols = gene_table[1].values

    # Create AnnData object
    var = pd.DataFrame(index=gene_symbols)
    adata = sc.AnnData(X=matrix, var=var, obs=pd.DataFrame(index=barcodes))
    print(f"Adata from {sample} created: {adata.shape}")

    # Save duplicate gene names
    duplicated_mask = pd.Series(gene_symbols).duplicated(keep=False)
    duplicate_genes = pd.Series(gene_symbols)[duplicated_mask]
    duplicates_file = os.path.join(data_dir, f"duplicated_genes_{cohort}_{sample}.csv")
    duplicate_genes.to_csv(duplicates_file, index=False, header=["duplicate_gene"])
    print(f"Found {duplicate_genes.nunique()} duplicated genes. Saved to {duplicates_file}")

    # Make gene names unique
    adata.var_names_make_unique()
    duplicates_after = adata.var_names[adata.var_names.duplicated()]
    if len(duplicates_after) == 0:
        print("All gene names are now unique!")
    else:
        print("There are still duplicates:", duplicates_after)

    # Add QC and metadata
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata.obs['case_id'] = sample
    adata.obs['Cohort'] = cohort
    #adata.obs['Diagnosis'] = diagnosis
    #adata.obs['IDH'] = idh
    #adata.obs['Primary_Recurrent'] =pt
    adata.obs['cell_id'] = barcodes
    adata.obs_names = [f"{cell}{study_id}_{sample}" for cell in adata.obs_names]


    # Save AnnData
    adata.write(output_file)
    print(f"Saved {sample} AnnData to {output_file}")
