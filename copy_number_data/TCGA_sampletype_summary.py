# ============================================================
# Script: TCGA CNA Sample Table Summary
# ============================================================
#
# Overview:
# This script performs summarization of TCGA
# copy number alteration (CNA) sample metadata, focusing on
# sample type distribution.
#
# ------------------------------------------------------------
# Objectives:
# ------------------------------------------------------------
# 1. Summarize distribution of TCGA sample types
# 2. Evaluate sample types at patient level
# 3. Identify patients with multiple primary tumor samples
# 4. Export patient-level summary table
#
# ------------------------------------------------------------
# Input:
# ------------------------------------------------------------
# - TCGA sample metadata table (CNA_Sampletable.tsv)
#   containing:
#   - case/sample identifiers
#   - sample type annotations (e.g., Primary Tumor, Normal Tissue)
#
# ------------------------------------------------------------
# Workflow:
# ------------------------------------------------------------
# 1. Load TCGA sample annotation table
# 2. Compute global distribution of sample types
# 3. Extract patient ID from TCGA barcode (first 12 characters)
# 4. Summarize unique sample types per patient
# 5. Identify patients with multiple primary tumor samples
# 6. Export patient-level summary table
#
# ------------------------------------------------------------
# Outputs:
# ------------------------------------------------------------
# - Console summary of sample type distribution
# - Patient-level sample type summary table
# - List of patients with >1 primary tumor sample
#
# ============================================================

import pandas as pd

# ------------------------------------------------------------
# Load TCGA CNA sample metadata table
# ------------------------------------------------------------
# This file contains TCGA case/sample-level annotations including:
# - sample type (e.g., Primary Tumor, Solid Tissue Normal)
# - case IDs (patient-level identifiers embedded in TCGA barcode)

CNA_sampletable = pd.read_csv( "/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/sampletable/CNA_Sampletable.tsv", sep = "\t") 

# ------------------------------------------------------------
# Overview of sample type distribution
# ------------------------------------------------------------
# This prints how many samples exist per sample category
# (e.g., Primary Tumor, Normal Tissue, etc.)
print(CNA_sampletable['sample_type'].value_counts())

# ------------------------------------------------------------
# Patient-level aggregation of sample types
# ------------------------------------------------------------
# TCGA sample barcodes contain patient identifiers in the first 12 characters.
# We extract this to group samples at patient level.
# Create a patient_id column from the first 12 characters of 'cases' 
CNA_sampletable['patient_id'] = CNA_sampletable['cases'].str[:12] 
 
# Group by patient_id and get the unique sample types per patient 
patient_summary = CNA_sampletable.groupby('patient_id')['sample_type'].unique().reset_index() 

# View the first few rows 
print(patient_summary.head())  

# Save patient-level summary for downstream analysis
# (useful for checking sample redundancy and dataset structure)
patient_summary.to_csv( 
    "/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/sampletable/CNV_patient_sample_summary.tsv", 
    sep="\t", index=False) 

# ------------------------------------------------------------
# Identify patients with multiple primary tumor samples
# ------------------------------------------------------------

# Count how many primary tumor samples each patient has
primary_counts = CNA_sampletable[CNA_sampletable['sample_type'] == 'Primary Tumor'] \
                 .groupby(CNA_sampletable['cases'].str[:12]) \
                 .size() \
                 .reset_index(name='primary_tumor_count')

# Identify patients with more than one primary tumor sample
multiple_primary = primary_counts[primary_counts['primary_tumor_count'] > 1]

print(multiple_primary)