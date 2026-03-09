import os
import csv
import itertools
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad('/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/phenotyping/adata_final.h5ad')
markers ='/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/annotations_markers/markers_poster.csv'
dataset_char = 'markers_adata_final_v3_marker_poster'
layer="log_norm"

file = open(markers, "r")
markers_list = list(csv.reader(file, delimiter=","))
file.close()

# Flatten the list using itertools.chain
markers_list = list(itertools.chain.from_iterable(markers_list))

# Find markers that are in adata.var
print('Finding marker genes in adata.var...')
markers_in_adata = [gene for gene in markers_list if gene in adata.var_names]
print(len(markers_in_adata))

# Set working directory (to save dotplots in desired folder)
print("Setting working directory...")
os.chdir("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/Dotplot")

# List different ways of groupby
groups = ['lv1_celltype']#['Cohort', 'cell_type', 'cell_type2', 'leiden', 'level0_celltype','level1_celltype']

#Create dotplots
print("Creating dotplots...")
for group in groups:
  if group in adata.var_names or group in adata.obs.columns:
      sc.set_figure_params(scanpy=True, fontsize=14)
      sc.pl.dotplot(adata, markers_in_adata, show=False, groupby=group, layer='log_norm', save=f'_{dataset_char}_{group}.png', use_raw=False)
  else:
    print(f"Warning: {group} not found in adata.var_names or adata.obs.columns.")
