"""
Script: Generate UMAP Visualizations

Description:
This script generates UMAP visualizations of an annotated single-cell RNA-seq
dataset using Scanpy. UMAPs can be colored by cell-type annotations,
sample metadata, clustering results, or gene expression.

The script:
1. Loads the annotated AnnData object.
2. Defines the annotation(s) or feature(s) to visualize.
3. Applies a consistent color palette for annotated cell types.
4. Generates UMAP plots for each selected annotation.
5. Saves each figure in both PNG and SVG formats.

The script can easily be adapted by modifying:
- the annotation(s) to plot,
- the color palette,
- the output directory,
- the AnnData object.

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""
import os
import scanpy as sc
import matplotlib.pyplot as plt

print('This is running')

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

# Load annotated single-cell dataset
print("Loading adata...")
adata = sc.read_h5ad("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/final_clean_sc_atlas/adata_final_cleaned.h5ad")

# Label used in output filenames
dataset_char = 'final_clean'

#-------------------------------------------------------------------------------
# 2. Set Output Directory
#-------------------------------------------------------------------------------

print("Setting working directory...")
os.chdir("/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/plots/UMAP/clustering/BBKNN/Leiden_1/FINAL/v3_thesis")

#-------------------------------------------------------------------------------
# 3. Select Annotations to Visualize
#-------------------------------------------------------------------------------

# Examples:
# attributes = ['SD3D', 'CD8A', 'FOXP3', 'IL7R']
# attributes = ['level1_celltype', 'level0_celltype', 'Diagnosis_label', 'Cohort', 'leiden']
# attributes = ['lv1_celltype', 'lv2_celltype', 'Cohort', 'kept_labels', 'batch']
# attributes = ['batch'

# Annotation(s) to display on the UMAP
attributes = ['level2_celltype', 'level3_celltype']

#-------------------------------------------------------------------------------
# 4. Define Cell-Type Color Palette
#-------------------------------------------------------------------------------

## Define a color dictionary
color_map = {'B_cells':'royalblue',
              'CD4_Tcell':'darkorange',
              'CD8_Tcell':'mediumseagreen',
              'Endothelial':'crimson',
              'Fibroblast':'darkviolet',
              'Malignant':'sienna',
              'Monocyte':'yellowgreen',
              'NK_cell':'mediumvioletred',
              'Oligodendrocyte':'mediumturquoise',
              'Pericyte':'aqua',
              'Plasma_B':'gray',
              'TAM':'lightpink',
              'T_reg':'mediumpurple',
              'cDC':'orangered',
              'pDC':'violet',
              'Microglia':'hotpink',
              'Macrophage':'lightskyblue',
              'Neutrophil' : 'gold'
              }

#-------------------------------------------------------------------------------
# 5. Generate UMAPs
#-------------------------------------------------------------------------------

print("Creating umaps...")

for attribute in attributes:
  # Check whether the requested annotation exists
  if attribute in adata.raw.var_names or attribute in adata.obs.columns:
    if attribute == 'leiden':
      sc.pl.umap(adata, color=attribute, palette = color_map, show=False, save=f'_{dataset_char}_{attribute}.png', legend_loc='on data', use_raw=True)
      sc.pl.umap(adata, color=attribute, palette = color_map, show=False, save=f'_{dataset_char}_{attribute}.svg', legend_loc='on data', use_raw=True)
    else:
      sc.pl.umap(adata, color=attribute, palette = color_map, show=False, save=f'_{dataset_char}_{attribute}.png', use_raw=True)
      sc.pl.umap(adata, color=attribute, palette = color_map, show=False, save=f'_{dataset_char}_{attribute}.svg', use_raw=True)
  else:
    print(f"Warning: {attribute} not found in adata.var_names or adata.obs.columns.")

#-------------------------------------------------------------------------------
# Notes
#-------------------------------------------------------------------------------

# Consider using use_raw=False when plotting gene expression if the desired
# values are stored in the log_norm layer rather than in adata.raw.

# If plotting metadata instead of cell-type annotations, remove the palette
# argument if a custom color map is not required.


# Example: plot expression of a single gene
#sc.pl.umap(adata, color = ['TMBS1X'])


# Example: create an output directory in python 
#import os
#os.makedirs("/net/beegfs/users/P086608/scRNA_glioma/Processing_scRNA/plots_practice/", exist_ok=True)
#sc.settings.figdir = "/net/beegfs/users/P086608/scRNA_glioma/Processing_scRNA/plots_practice/"
#sc.pl.umap(adata, color=['IL1R2'], save='_IL1R2.png', show=False)
