import os
import scanpy as sc
import matplotlib.pyplot as plt

print('This is running')

# -------------------------------------------------
# Load adata
# -------------------------------------------------
print("Loading adata...")
adata = sc.read_h5ad("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/phenotyping/adata_final.h5ad")

dataset_char = 'all_final'

# -------------------------------------------------
# Set figure directory
# -------------------------------------------------
print("Setting working directory...")
output_dir = "/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/UMAP/clustering/BBKNN/Leiden_1/FINAL"
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir



print("Creating Fibroblast highlight UMAP...")

celltype_column = 'lv1_celltype'  # adjust if needed
target_celltype = 'nan'        # must match category exactly

if celltype_column in adata.obs.columns:

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
