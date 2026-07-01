# =============================================================================
# Script: Add GEO-Derived Metadata into AnnData Object (GSE278456)
# =============================================================================
#
# Description:
# This script processes and integrates sample-level metadata from a GEO
# metadata text file into an existing single-cell RNA-seq AnnData object.
#
# It standardizes clinical annotations (e.g., diagnosis, IDH status) and
# merges them with per-cell metadata using the shared sample identifier
# (case_id). The final AnnData object contains both expression data and
# enriched biological annotations for downstream analysis.
#
# Workflow:
# 1. Define paths to GEO metadata file and AnnData object.
# 2. Extract relevant sample-level metadata fields from GEO file.
# 3. Convert raw metadata into a structured pandas DataFrame.
# 4. Parse and clean key clinical variables:
#    - Extract diagnosis from tumor classification
#    - Standardize IDH status (WT / Mut)
# 5. Load existing AnnData object.
# 6. Preserve original cell index before merging.
# 7. Merge metadata into adata.obs using `case_id`.
# 8. Restore original AnnData indexing.
# 9. Save updated AnnData object with integrated metadata.
#
# Inputs:
# - GEO metadata file (GSE278456_meta.txt)
# - Combined AnnData object (.h5ad)
#
# Outputs:
# - Updated AnnData object with additional `.obs` metadata columns
#
# Notes:
# - `case_id` is used as the primary join key between metadata and AnnData.
# - IDH status is simplified into WT / Mut for consistency across analyses.
# =============================================================================

import pandas as pd
import scanpy as sc

# Define file paths for metadata and AnnData object
meta_path = "/net/beegfs/users/P086608/scRNA_glioma/GSE278456/meta/GSE278456_meta.txt"
adata_path = "/net/beegfs/users/P086608/scRNA_glioma/GSE278456/adata_files/combined_raw2_GSE278456.h5ad"

# Define metadata fields to extract from GEO file and rename consistently
target_fields = {
    "!Sample_title": "case_id",
    "!Tumor_type": "tumor_type",
    #"!Tumor_grade": "tumor_grade",
    "!IDH_status": "IDH"
}

# Initialize dict to store extracted metadata columns
data = {v: None for v in target_fields.values()}

# Read metadata file
with open(meta_path, "r") as f:
    for line in f:
        for raw_key, clean_key in target_fields.items():
            if line.startswith(raw_key):
                values = line.strip().split("\t")[1:]          # skip the first column
                values = [v.replace('"', '') for v in values]  # remove quotes
                data[clean_key] = values

# Convert to DataFrame
df = pd.DataFrame(data)

# Extract diagnosis and IDH from tumor_type
df['Diagnosis'] = df['tumor_type'].str.extract(r'classification: (\w+)', expand=False)

# Clean up strings
# Remove old column
df.drop(columns=['tumor_type'], inplace=True)

# Clean IDH field by removing prefix text
df['IDH'] = df['IDH'].str.replace('idh.status: ', '', regex=False)

# Standardize IDH labels for downstream consistency
df['IDH'] = df['IDH'].map({
    'wild type': 'WT',
    'mutant': 'Mut'
})

# ------------------------------------------------------------------
# Load AnnData object containing single-cell expression matrix
# -----------------------------------------------------------------
adata = sc.read_h5ad(adata_path)

# Create column that will later become index
adata.obs['later_index'] = adata.obs.index.tolist()
print(adata)
print(adata.obs)

# -------------------------------------------------------------------
# Merge sample-level metadata into per-cell AnnData observations
# -----------------------------------------------------------------
adata.obs = adata.obs.merge(
    df,
    on="case_id",   # must match the sample column in adata.obs
    how="left"
)

# Now set 'later_index' as the index
adata.obs.set_index('later_index', inplace=True)
adata.obs.index.name = None


print(adata)
print(adata.obs)

# -----------------------------------------------------
# Save updated AnnData object with merged metadata
# -----------------------------------------------------
adata.write_h5ad(adata_path)

print("Metadata merged successfully! Updated adata.obs columns:")
print(adata.obs.columns)