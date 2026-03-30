"""
Script: Visualization of Pseudobulk Cell Type Fractions by Diagnosis

Description:
This script visualizes cell-type compositions derived from pseudobulk 
deconvolution of single-cell RNA-seq data, stratified by patient diagnosis.

The script:
1. Loads single-cell RNA-seq AnnData object containing sample and diagnosis metadata.
2. Loads pseudobulk deconvolution results (Statescope-derived cell fractions).
3. Extracts and maps diagnosis labels from the AnnData object.
4. Merges diagnosis metadata into the deconvolution output.
5. Converts the dataset from wide to long format for plotting with Seaborn.
6. Creates boxplots for each cell type showing distributions by diagnosis:
   - Boxplots summarize group distributions.
   - Stripplots overlay individual sample points for transparency.
7. Saves the figure for downstream reporting and analysis.

Author: Mabel Pronk
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

#-------------------------------------------------------------------------------
# 1. Load single-cell RNA-seq data (AnnData)
#-------------------------------------------------------------------------------
adata = sc.read_h5ad('/net/beegfs/users/P086608/scRNA_glioma/new_approach_processing_scRNA/data/final_clean_sc_atlas/adata_final_raw.h5ad')

#-------------------------------------------------------------------------------
# 2. Load Pseudobulk Deconvolution Results
#-------------------------------------------------------------------------------
deconv_path = '/net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_0001_rep_10/fractions3.csv'
df_deconv = pd.read_csv(deconv_path, index_col=0)

#-------------------------------------------------------------------------------
# 3. Extract Diagnosis Mapping from AnnData
#-------------------------------------------------------------------------------
# Pull sample ID and diagnosis label from AnnData observations
# Remove duplicates to get a single entry per sample
diagnosis_map = adata.obs[['sample_ID', 'Diagnosis_label']].set_index('sample_ID')

#-------------------------------------------------------------------------------
# 4. Merge Diagnosis Labels into Deconvolution Results
#-------------------------------------------------------------------------------
df_plot = df_deconv.join(diagnosis_map).dropna(subset=['Diagnosis_label'])

#-------------------------------------------------------------------------------
# 5. Melt Data from Wide to Long Format for Seaborn
#-------------------------------------------------------------------------------
# Converts columns [Sample | Malignant | T-cell | Diagnosis] into:
# [Sample | Diagnosis | Cell_Type | Fraction]
df_long = df_plot.reset_index().melt(
    id_vars=['index', 'Diagnosis_label'], 
    var_name='Cell_Type', 
    value_name='Fraction'
)

#-------------------------------------------------------------------------------
# 6. Create Boxplots
#-------------------------------------------------------------------------------
plt.figure(figsize=(16, 8))

# Boxplot per cell type stratified by diagnosis
sns.boxplot(
    data=df_long, 
    x='Cell_Type', 
    y='Fraction', 
    hue='Diagnosis_label', 
    palette='Set2',
    showfliers=False  # hides outliers for cleaner visualization
)

# Overlay individual sample points (jittered)
sns.stripplot(
    data=df_long, 
    x='Cell_Type', 
    y='Fraction', 
    hue='Diagnosis_label', 
    dodge=True, 
    color='black', 
    alpha=0.3, 
    size=2,
    legend=False
)

# Customize plot appearance
plt.title('Cell Type Fractions by Diagnosis (Pseudobulk)', fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Fraction')
plt.tight_layout()

#-------------------------------------------------------------------------------
# 7. Save Figure
#-------------------------------------------------------------------------------
plt.savefig(
    '/net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_0001_rep_10/boxplot_cell_types_by_diagnosis.png', 
    dpi=300
)
plt.show()
