import os
import scanpy as sc
import matplotlib.pyplot as plt

print('This is running')
# Load adata
print("Loading adata...")
adata = sc.read_h5ad("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/final_clean_sc_atlas/adata_final_cleaned.h5ad")
dataset_char = 'final_clean'

# Set working directory (to save umaps in desired folder)
print("Setting working directory...")
os.chdir("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/UMAP/clustering/BBKNN/Leiden_1/FINAL/v3_cleaned")


#attributes = ['SD3D','CD8A', 'FOXP3', 'IL7R']
#attributes = ['level1_celltype', 'level0_celltype', 'Diagnosis_label', 'Cohort','leiden']
#attributes = ['lv1_celltype','lv2_celltype', 'Cohort', 'kept_labels', 'batch']
#attributes = ['batch']
attributes = ['level2_celltype', 'level3_celltype', 'Cohort']
#Create umaps
print("Creating umaps...")
for attribute in attributes:
  if attribute in adata.raw.var_names or attribute in adata.obs.columns:
    if attribute == 'leiden':
      sc.pl.umap(adata, color=attribute, show=False, save=f'_{dataset_char}_{attribute}.png', legend_loc='on data', use_raw=True)
    else:
      sc.pl.umap(adata, color=attribute, show=False, save=f'_{dataset_char}_{attribute}.png', use_raw=True)
  else:
    print(f"Warning: {attribute} not found in adata.var_names or adata.obs.columns.")

# Should probably use use_raw = FALSE as we want to plot log_norm values which are in the log_norm layer 
