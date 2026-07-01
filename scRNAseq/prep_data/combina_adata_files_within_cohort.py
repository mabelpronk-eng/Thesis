# ============================================================
# Script: Merge Individual AnnData Files into a Combined Dataset
# ============================================================
#
# Description:
# This script combines multiple sample-level AnnData (.h5ad) files from a
# single scRNA-seq cohort into one merged AnnData object. Prior to merging,
# the script identifies the set of genes shared across all samples to ensure
# a consistent gene expression matrix.
#
# Workflow:
# 1. Identify all sample-level AnnData files in the input directory.
# 2. Load each AnnData object into memory.
# 3. Determine the intersection of genes present across all samples.
# 4. Subset each AnnData object to the shared gene set.
# 5. Concatenate all samples into a single AnnData object.
# 6. Inspect the merged dataset.
# 7. Save the combined AnnData object for downstream analyses.
#
# Inputs:
# - Individual sample-level .h5ad files
#
# Outputs:
# - Combined AnnData (.h5ad) containing all cells from the cohort
#
# Purpose:
# To generate a unified scRNA-seq dataset from multiple processed samples,
# providing a standardized input for downstream quality control, normalization,
# integration, clustering, and cell type annotation.
# ============================================================
import os
import anndata as ad
import scanpy as sc


folder_path = "/net/beegfs/users/P086608/scRNA_glioma/GSE222522/adata_files"   # Path to the folder containing the files
file_starting = 'GSE222522_'  # File ending to recognize relevant files in folder
output_file =  '/net/beegfs/users/P086608/scRNA_glioma/GSE222522/adata_files/combined_raw2_GSE222522.h5ad' # Directory and name of the output file

#-----------------------------------------------------------------------------
# 1. Get list of all the samples
#-----------------------------------------------------------------------------

# List all files starting with file_starting
file_list = [f for f in os.listdir(folder_path) if f.startswith(file_starting)]
#print(file_list)
print(len(file_list))

#-----------------------------------------------------------------------------
# 2. Load all AnnData objects
#-----------------------------------------------------------------------------
adatas = [sc.read_h5ad(os.path.join(folder_path, f)) for f in file_list]

#-----------------------------------------------------------------------------
# 3. Ensure common genes across all AnnData objects
#-----------------------------------------------------------------------------
print('Finding common genes...')
common_genes = set(adatas[0].var_names)
for adata in adatas[1:]:
    common_genes &= set(adata.var_names)
common_genes = sorted(list(common_genes))  # Sort for consistency
print('Number of common genes:', len(common_genes))

# Subset each AnnData to only include the common genes
adatas = [adata[:, common_genes] for adata in adatas]

#-----------------------------------------------------------------------------
# 4. Combine all AnnData objects
#-----------------------------------------------------------------------------
print('Combining adata...')
combined_adata = ad.concat(
    adatas,
    axis=0,           # CONCATENATE CELLS (not genes)
    join='inner',     # Only keep shared genes (redundant since we already subset)
    index_unique=None # Keep original obs_names (make unique manually if needed)
)

#-----------------------------------------------------------------------------
# 5. Inspect adata file
#-----------------------------------------------------------------------------
# Inspect the combined AnnData
print(combined_adata)
print(combined_adata.var)
print(combined_adata.obs)

# Save the combined AnnData
combined_adata.write(output_file) 