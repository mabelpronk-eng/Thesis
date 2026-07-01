"""
Script: Generate Highlighted UMAP for a Selected Cell Population

Description:
This script generates a UMAP highlighting a single selected cell population
within an annotated single-cell RNA-seq dataset. The selected population is
displayed in a specified color, while all remaining cells are shown in grey,
allowing clear visualization of its spatial distribution.

The script:
1. Loads the annotated AnnData object.
2. Creates the output directory if it does not exist.
3. Selects a cell annotation and target cell population.
4. Generates a UMAP highlighting only the selected population.
5. Saves the resulting figure.

The script can easily be adapted by modifying:
- the annotation column,
- the target cell population,
- the highlight color,
- the output directory.
"""
import os
import scanpy as sc
import matplotlib.pyplot as plt

print('This is running')

# -------------------------------------------------
# Load Data
# -------------------------------------------------
# Load annotated single-cell dataset
print("Loading adata...")
adata = sc.read_h5ad("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/phenotyping/adata_final.h5ad")

# Label used in output filename
dataset_char = 'all_final'

#-------------------------------------------------------------------------------
# 2. Set Output Directory
#-------------------------------------------------------------------------------
print("Setting working directory...")
output_dir = "/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/UMAP/clustering/BBKNN/Leiden_1/FINAL"
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

#-------------------------------------------------------------------------------
# 3. Select Cell Population to Highlight
#-------------------------------------------------------------------------------

print(f"Creating {target_celltype} highlight UMAP...")

# Cell annotation containing the desired labels
celltype_column = 'lv1_celltype'  # adjust if needed

# Cell population to highlight
target_celltype = 'nan'        # must match category exactly

#-------------------------------------------------------------------------------
# 4. Generate Highlighted UMAP
#-------------------------------------------------------------------------------
# Check whether the annotation exists
if celltype_column in adata.obs.columns:

    # Check whether the selected cell population exists
    if target_celltype in adata.obs[celltype_column].unique():

        sc.pl.umap(
            adata,
            color=celltype_column,
            groups=[target_celltype],     # highlight only this group
            palette=['red'],              # color for Fibroblast
            na_color='lightgrey',         # all other cells grey
            show=False,
            save=f'_{dataset_char}_{target_celltype}_highlight.png'
        )

    else:
        print(f"{target_celltype} not found in {celltype_column}")

else:
    print(f"{celltype_column} not found in adata.obs")


print("Done.")
