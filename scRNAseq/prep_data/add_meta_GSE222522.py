# =============================================================================
# Script: Add GEO Sample Metadata to AnnData Object (GSE222522)
# =============================================================================
#
# Description:
# This script extracts clinical and sample-level metadata from a GEO
# series matrix file and merges it into an existing AnnData object.
#
# The metadata includes diagnostic labels and molecular characteristics
# (e.g., IDH mutation status, primary vs recurrent tumor status),
# which are linked to single-cell samples via GSM identifiers.
#
# Workflow:
# 1. Load GEO series matrix file and extract sample (GSM) identifiers.
# 2. Parse sample-level clinical annotations from GEO characteristics fields.
# 3. Standardize key clinical variables (e.g., IDH status).
# 4. Load existing AnnData object.
# 5. Extract GSM identifiers from AnnData sample metadata.
# 6. Merge GEO metadata with AnnData observations using GSM as key.
# 7. Restore original cell indexing structure.
# 8. Save updated AnnData object with integrated metadata.
#
# Inputs:
# - GEO series matrix file (GSE222522-GPL24676_series_matrix.txt)
# - Existing AnnData file (.h5ad)
#
# Outputs:
# - Updated AnnData object with added metadata columns in `.obs`
#
#
# Notes:
# - GSM IDs are used as the linking key between GEO metadata and AnnData.
# - IDH status is harmonized into simplified categories (WT / Mut).
#
# Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# =============================================================================

import pandas as pd
import re
import scanpy as sc

adata_path = '/net/beegfs/users/P086608/scRNA_glioma/GSE222522/adata_files/combined_raw2_GSE222522.h5ad'

# -------------------------------------------------------------------------
# Load GEO series matrix file containing sample-level clinical metadata
# -------------------------------------------------------------------------
with open("/net/beegfs/users/P086608/scRNA_glioma/GSE222522/meta/GSE222522-GPL24676_series_matrix (1).txt") as f:
    lines = f.readlines()

# -------------------------------------------------------------------------
# Extract GSM identifiers (GEO sample IDs)
# -------------------------------------------------------------------------
for line in lines:
    if line.startswith("!Sample_geo_accession"):
        gsm_list = re.findall(r'GSM\d+', line)
        break

# Create metadata table starting from GSM identifiers
meta = pd.DataFrame({"GSM": gsm_list})

# -------------------------------------------------------------------------
# Function to extract clinical/sample characteristics from GEO file
# -------------------------------------------------------------------------
def extract_characteristic(lines, key):
    for line in lines:
        # Look for lines containing sample characteristics
        if line.startswith("!Sample_characteristics_ch1") and key in line:
            # Extract quoted values from GEO metadata line
            values = re.findall(r'"([^"]+)"', line)
            # Split key:value format and keep only the value part
            return [v.split(": ", 1)[1] for v in values]
    return None

# Extract clinical annotations from GEO metadata
meta["Diagnosis_label"] = extract_characteristic(lines, "diagnosis:")
meta["IDH"] = extract_characteristic(lines, "genotype:")
meta["Primary_Recurrent"] = extract_characteristic(lines, "primary/recurrent:")

# -------------------------------------------------------------------------
# Standardize IDH labels for downstream analysis consistency
# -------------------------------------------------------------------------
meta["IDH"] = meta["IDH"].replace({
    "IDH Wild Type": "WT",
    "IDH Mutant": "Mut",
})
meta.head()

# -------------------------------------------------------------------------
# Load existing AnnData object
# -------------------------------------------------------------------------
adata = sc.read_h5ad(adata_path)

# -------------------------------------------------------------------------
# Preserve original cell index temporarily (needed after merge)
# -------------------------------------------------------------------------
adata.obs['later_index'] = adata.obs.index.tolist()
print(adata)
print(adata.obs)

# Extract GSM IDs from AnnData sample metadata (case-level annotation)
adata.obs["GSM"] = adata.obs["case_id"].str.extract(r"(GSM\d+)")

# Merge GEO metadata into AnnData observations using GSM as key
adata.obs = adata.obs.merge(
    meta,
    on="GSM",
    how="left"
)

# Now set 'later_index' as the index
# Restore original cell indexing after merge
adata.obs.set_index('later_index', inplace=True)
adata.obs.index.name = None


print(adata)
print(adata.obs)

# -------------------------------------------------------------------------
# Save updated AnnData object with merged metadata
# -------------------------------------------------------------------------
adata.write_h5ad(adata_path)

print("Metadata merged successfully! Updated adata.obs columns:")
print(adata.obs.columns)
