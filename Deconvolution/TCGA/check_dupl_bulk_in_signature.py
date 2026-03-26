"""
Script: Identification of Duplicate Genes in Bulk RNA-seq and Cross-Reference with Signature Markers

Description:
This script identifies duplicated gene entries in a TCGA bulk RNA-seq dataset and 
cross-references them with marker genes from a predefined cell type signature matrix. 
It reports any duplicated genes that are also annotated as marker genes, as these may 
affect downstream deconvolution analyses. Additionally, it provides a preview of the 
expression values for selected problematic genes to assess potential discrepancies.

Author: Mabel Pronk
"""

import pandas as pd

# 1. Load input files
# Bulk RNA-seq expression matrix (genes x samples)
bulk_path = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/input/transcriptome_matrix_TCGA_clean_column.csv"

# Signature matrix containing marker gene annotations
sig_path = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/input/signature_17celltypes_GBM.txt"

bulk = pd.read_csv(bulk_path, index_col=0)
sig = pd.read_csv(sig_path, index_col=0, sep="\t")

# 2. Identify duplicated gene entries in the bulk dataset
# keep=False marks all occurrences of duplicated gene names
is_duplicate = bulk.index.duplicated(keep=False)
duplicate_names = bulk.index[is_duplicate].unique()

# 3. Cross-reference duplicates with marker genes from the signature matrix
# Identify genes that are both duplicated and annotated as markers (IsMarker == True)
marker_genes = sig[sig['IsMarker'] == True].index
overlap_problem = [gene for gene in duplicate_names if gene in marker_genes]

# 4. Report results and inspect potential issues
if overlap_problem:
    print(f"Found {len(overlap_problem)} duplicate bulk genes that are MARKERS in your signature.")
    print("--- Problematic Marker Genes ---")
    print(overlap_problem)
    
    # Preview expression values for the first few problematic genes
    # This helps assess whether duplicate rows contain inconsistent values
    print("\nPreview of bulk values for these markers (checking for variance):")
    for gene in overlap_problem[:3]:
        print(f"\nGene: {gene}")
        print(bulk.loc[[gene]].iloc[:, :5])  # Display first 5 samples
else:
    print("Clean! None of your bulk duplicates are listed as 'IsMarker' in the signature.")
