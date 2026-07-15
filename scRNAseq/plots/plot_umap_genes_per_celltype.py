"""
Script: Generate UMAPs for Marker Gene Expression

Description:
This script generates UMAP visualizations showing the expression of selected
marker genes within an annotated single-cell RNA-seq dataset. Marker genes are
organized by cell type, and each gene is plotted individually. Figures are
saved into separate folders corresponding to each cell type.

The script:
1. Loads the annotated AnnData object.
2. Defines marker genes for one or more cell populations.
3. Creates an output directory for each cell type.
4. Checks whether each marker gene is present in the dataset.
5. Generates a UMAP displaying normalized expression of each marker gene.
6. Saves all figures to their respective output folders.

The script can easily be adapted by modifying:
- the marker gene dictionary,
- the expression layer,
- the output directory,
- the AnnData object.

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""
import os
import scanpy as sc
import matplotlib.pyplot as plt

print("This is running")

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

# Load annotated single-cell dataset
print("Loading adata...")
adata = sc.read_h5ad("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/phenotyping/adata_final_with_sampleID.h5ad")

# Label used in output filenames
dataset_char = 'final'

#-------------------------------------------------------------------------------
# 2. Set Output Directory
#-------------------------------------------------------------------------------

# Define directory in which all figures will be saved
print("Setting working directory...")
base_dir = "/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/UMAP/clustering/BBKNN/Leiden_1/FINAL/v3"
os.chdir(base_dir)

#-------------------------------------------------------------------------------
# 3. Define Marker Genes
#-------------------------------------------------------------------------------

# Marker genes grouped by cell type
marker_genes = {
    #"B_cells": ["CD79A", "MS4A1"],
    #"CD4_Tcell": ["CD3D",'CD3G', "IL7R", 'CD4'],
    #"CD8_Tcell": ["CD3D", "GZMK", "CD8A", 'CD8B'],
    #"Endothelial": ["CLDN5", "VWF", 'CD34'],
    #"Fibroblast": ["COL1A1", "FBLN2", "DCN", "VCAN"],
    #"Macrophage": ["C1QA", "FCGR1A", "CD68", "CD14", 'CD163', 'CLEC4E'],
    "Malignant": ["SOX2", 'EGFR','MKI67','TP53INP1'],
    #"Microglia": ["P2RY12", "TMEM119", 'CX3CR1'],
    #"Monocyte": ["CD14", "VCAN", "LYZ", "S100A8", 'FCN1'],
    #"NK_cell": ["GNLY", "NKG7", "KLRC1", 'NCAM1'],
    #"Oligodendrocyte": ["MAG", "MOG", "PLP1", "MBP"],
    #"Pericyte": ["RGS5", "PDGFRB", "DCN"],  # corrected gene PDGFRB
    #"Plasma_B": ["CD79A", 'CD38','MZB1','IRF4','SLAMF7'],
    #"T_reg": ["CD3D", "FOXP3", "CTLA4", "IL2RA"],
    #"cDC": ["FCER1A", "CD1C",'HLA-DRA', 'CLEC9A', 'ITGAX','CLEC12A'],
    #'matureDC':['CCR7', 'LAMP3', 'CD83'], 
    #"pDC": ["LILRA4", "TPM4", 'TPM2','GZMB', 'IRF4'],
    #"Neutrophils": ["IL1R2", "CXCR2", "FPR2", 'FCGR3B', 'ITGAM', 'SELL']
}

# Alternative marker gene panel (currently unused)
marker_dict = {
    # Classical dendritic cells
    "cDC1": ["CLEC9A", "CADM1", "IDO1", "CST3"],
    "cDC2": ["CD1C", "CLEC10A", "FCER1A", "HLA-DRA"],
    "matureDC": ["LAMP3", "CCR7", "IL4I1", "CD83"],
    "pDC": ["LILRA4", "IL3RA", "GZMB", "TCF4"],

    # Monocyte / macrophage states
    "Monocyte": ["CD14", "FCN1", "S100A9", "VCAN"],
    "TAM": ["CD163", "MAFB", "CCL3", "TREM2"],
    "Microglia": ["C1QA", "C1QB", "C1QC", "CX3CR1", "APOE"],

    # Lymphoid contamination
    "T_cell": ["CD3D", "CD3E", "NKG7", "CCL5"]
}


#-------------------------------------------------------------------------------
# 4. Generate Marker Gene UMAPs
#-------------------------------------------------------------------------------

print("Creating UMAPs for marker genes...")

for celltype, genes in marker_genes.items():
    
    # Make folder for this cell type
    folder = os.path.join(base_dir, celltype)
    os.makedirs(folder, exist_ok=True)
    # VERY IMPORTANT: set Scanpy figdir to absolute path
    sc.settings.figdir = folder

    print(f"\nProcessing {celltype} → saving to {folder}")

    # Loop through genes
    for gene in genes:

        # Check if gene exists in raw.var or var_names
        if gene not in adata.var_names:
            print(f"  ⚠️ Warning: {gene} not found in adata.var_names.")
            continue

        print(f"  Plotting {gene}...")

        # Save inside the celltype folder
        save_path = f"{dataset_char}_{celltype}_{gene}.png"

        sc.pl.umap(
            adata,
            color=gene,
            show=False,
            save=save_path,
            layer="log_norm"

        )

print("\nAll done.")
