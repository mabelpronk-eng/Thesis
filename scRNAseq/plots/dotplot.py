"""
Script: Generate Marker Gene Dotplots

Description:
This script generates Scanpy dotplots to visualize the expression of selected
marker genes across annotated cell populations in a single-cell RNA-seq dataset.

The script:
1. Loads the annotated AnnData object.
2. Reads a list of marker genes from a CSV file.
3. Retains only marker genes present in the dataset.
4. Iterates over one or more cell annotations (groupings).
5. Creates and saves a dotplot showing:
   - Dot size: fraction of cells expressing the gene.
   - Dot color: average normalized expression level.
6. Saves the resulting figures to the specified output directory.

The script can easily be adapted by modifying:
- the marker gene list,
- the grouping variable(s),
- the expression layer,
- the output directory.
"""

import os
import csv
import itertools
import scanpy as sc
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

# Load annotated single-cell dataset
adata = sc.read_h5ad('/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/final_clean_sc_atlas/adata_final_cleaned.h5ad')

# Marker gene list
markers ='/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/annotations_markers/markers_used_thesis.csv'

# Label used in output filenames
dataset_char = 'markers_refinement_thesis'

# Expression layer used for visualization
layer="log_norm"

#-------------------------------------------------------------------------------
# 2. Load Marker Genes
#-------------------------------------------------------------------------------

file = open(markers, "r")
markers_list = list(csv.reader(file, delimiter=","))
file.close()

# Flatten the list using itertools.chain
markers_list = list(itertools.chain.from_iterable(markers_list))

# Find markers that are in adata.var
print('Finding marker genes in adata.var...')
markers_in_adata = [gene for gene in markers_list if gene in adata.var_names]
print(len(markers_in_adata))

#-------------------------------------------------------------------------------
# 3. Set Output Directory
#-------------------------------------------------------------------------------

# Set working directory (to save dotplots in desired folder)
print("Setting working directory...")
os.chdir("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/Dotplot")

#-------------------------------------------------------------------------------
# 4. Generate Dotplots
#-------------------------------------------------------------------------------

# Annotation(s) used to group cells in the dotplot
groups = ['level3_celltype']
# Example alternatives:
# ['Cohort', 'cell_type', 'cell_type2', 'leiden',
#  'level0_celltype', 'level1_celltype']

#Create dotplots
print("Creating dotplots...")
for group in groups:
  if group in adata.var_names or group in adata.obs.columns:
      sc.set_figure_params(scanpy=True, fontsize=14)
      sc.pl.dotplot(adata, markers_in_adata, show=False, groupby=group, layer='log_norm', save=f'_{dataset_char}_{group}.png', use_raw=False)
      sc.pl.dotplot(adata, markers_in_adata, show=False, groupby=group, layer='log_norm', save=f'_{dataset_char}_{group}.pdf', use_raw=False)
  else:
    print(f"Warning: {group} not found in adata.var_names or adata.obs.columns.")
