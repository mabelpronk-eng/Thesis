# ============================================================
# Script: Convert Individual 10x Genomics HDF5 Files to AnnData
# ============================================================
#
# Description:
# This script processes multiple single-cell RNA-seq datasets stored as
# 10x Genomics HDF5 files and converts each sample into an individual
# AnnData (.h5ad) object. During processing, basic quality control
# metrics and sample metadata are added to each dataset.
#
# Workflow:
# 1. Identify all 10x HDF5 files in the input directory.
# 2. Load each HDF5 file into an AnnData object.
# 3. Ensure gene names are unique.
# 4. Calculate basic quality control metrics (e.g., mitochondrial content).
# 5. Add sample- and cohort-level metadata.
# 6. Generate unique cell identifiers by appending study and sample IDs.
# 7. Save each processed sample as an individual AnnData (.h5ad) file.
#
# Inputs:
# - Multiple 10x Genomics HDF5 files (*_filtered_counts.h5)
#
# Outputs:
# - One processed AnnData (.h5ad) file per sample
#
#
# Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# ============================================================

import os
import scanpy as sc

#-----------------------------------------------------------------------------
# 0. Prepare variables
#-----------------------------------------------------------------------------
folder_path = "/net/beegfs/users/P086608/scRNA_glioma/GSE278456"
output_dir  = "/net/beegfs/users/P086608/scRNA_glioma/GSE278456/adata_files"
os.makedirs(output_dir, exist_ok=True)
study_id    = "_GSE278456"
cohort      = "GSE278456"

# Get all .h5 files matching the naming pattern
file_list = [f for f in os.listdir(folder_path) if f.endswith("_filtered_counts.h5")]

print(f"Found {len(file_list)} files:")
for f in file_list:
    print(f"  {f}")

#-----------------------------------------------------------------------------
# 1. Process each file
#-----------------------------------------------------------------------------
for file_name in file_list:
    file_path = os.path.join(folder_path, file_name)

    # Parse identifiers from filename
    # Example: GSM8546881_GBM006_Myl_filtered_counts.h5
    parts = file_name.split('_')
    sample_geo_accession = parts[0]             # GSM8546881
    sample_name = parts[1]                      # GBM006
    filter = parts[2]                           # Myl 

    print(f"\nProcessing {file_name}")
    #print(f"Sample GEO accession: {sample_geo_accession}")
    #print(f"Sample name: {sample_name}")

    # Read data
    adata = sc.read_10x_h5(file_path)

    # Make gene names unique
    adata.var_names_make_unique()

    #-----------------------------------------------------------------------------
    # 2. Add metadata
    #-----------------------------------------------------------------------------
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Add sample and cohort info to obs
    adata.obs['case_id'] = sample_name
    adata.obs['geo_accession'] = sample_geo_accession
    adata.obs['Cohort'] = cohort
    adata.obs['filter'] = filter
    adata.obs['cell_id'] = adata.obs_names.copy()

    # Add study ID to obs names for uniqueness
    adata.obs_names = (adata.obs_names + study_id + '_' + sample_name)

    #-----------------------------------------------------------------------------
    # 3. Save processed AnnData
    #-----------------------------------------------------------------------------
    output_file = os.path.join(output_dir, f"{sample_geo_accession}_{sample_name}_{filter}.h5ad")
    adata.write(output_file)
    print(f"Saved to {output_file}")
