# ============================================================
# Script: Download and Extract GEO Single-Cell RNA-seq Dataset
# ============================================================
#
# Description:
# This script downloads a raw single-cell RNA-seq dataset from the
# NCBI GEO repository (GSE222522), extracts the compressed archive,
# and prepares the raw 10x Genomics files for downstream processing.
#
# The dataset is provided as a gzip-compressed ZIP archive (.zip.gz),
# which requires a two-step extraction process.
#
# Workflow:
# 1. Define working directory and ensure it exists.
# 2. Change into the working directory.
# 3. Download GEO dataset archive using wget.
# 4. Decompress the gzip layer (.gz → .zip).
# 5. Extract the ZIP archive to obtain raw files.
# 6. List extracted files for verification.
#
# Inputs:
# - GEO accession: GSE222522
# - Remote archive: .zip.gz file from NCBI GEO
#
# Outputs:
# - Extracted raw single-cell RNA-seq files in 10x format
#   (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz, etc.)
#
# Purpose:
# To retrieve and unpack raw GEO single-cell RNA-seq data so it can
# be processed into AnnData objects for downstream analysis.
# ============================================================

import os

workspace_folder = '/net/beegfs/users/P086608/scRNA_glioma/GSE222522'
os.makedirs(workspace_folder, exist_ok=True) #Make sure the folder exists and creates it if it doesn't 
os.chdir(workspace_folder) #changes directory to this 

### Download zip.gz ###
# Download URL and file names
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE222522&format=file"

# File names for intermediate and final extraction steps
gz_file = "GSE222522_raw_counts_matrix.zip.gz"
zip_file = "GSE222522_raw_counts_matrix.zip"

# ------------------------------------------------------------
# Step 1: Download compressed GEO archive
# ------------------------------------------------------------
os.system(f"wget -O {gz_file} '{url}'")

# ------------------------------------------------------------
# Step 2: Decompress gzip layer (.gz → .zip)
# ------------------------------------------------------------
os.system(f"gunzip -f {gz_file}")  # This will produce a .zip file

# ------------------------------------------------------------
# Step 3: Extract ZIP archive
# ------------------------------------------------------------
os.system(f"unzip -o {zip_file}")

print(f"✅ Download and extraction complete.\nFiles are in {workspace_folder}")

# List extracted files for confirmation
os.system("ls -lh")

### Download TAR (alternative GEO format) ###
# The following block is an alternative workflow for datasets
# distributed as TAR archives (e.g., GSE278456).
#url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE278456&format=file" #download link
#output_file = "GSE278456_RAW.tar" #name of downloaded file

#os.system(f"wget -O {output_file} '{url}'")
#os.system(f"tar -xvf {output_file}") #Extract all the data probably not handy to have put this in! First view data 

#print(f"Download and extraction complete. Files are in {workspace_folder}")
