"""
Script: Visualization of Deconvoluted Cell Type Fractions (Stacked Bar Plot)

Description:
This script visualizes cell type composition inferred by Statescope (or another deconvolution method) from bulk RNA-seq data.
Samples are sorted by malignant cell fraction in descending order to highlight gradients 
in tumor purity and microenvironment composition.

A stacked bar plot is generated, where each bar represents a sample and each segment 
corresponds to the relative fraction of a specific cell type. The plot is automatically 
scaled based on the number of samples and saved for downstream analysis and reporting.

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

#-------------------------------------------------------------------------------
# 1. Load the Data
#-------------------------------------------------------------------------------

# Path to Statescope deconvolution output (cell fractions per sample)
deconv_path = '/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output/fractions3.csv'
#deconv_path = '/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output/cibersort_fractions_TCGAbulk_gbm.tsv'

# Output directory and file name for the figure
save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output"
#save_dir = '/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output/plots'
save_name = "stacked_cell_fractions_sorted_by_malignant.png"

# Load data (rows = samples, columns = cell types)
df = pd.read_csv(deconv_path, index_col=0)

#-------------------------------------------------------------------------------
# 2. Sort and Generate Stacked Bar Plot
#-------------------------------------------------------------------------------

# Sort samples by malignant fraction (descending) to visualize purity gradient
df_sorted = df.sort_values(by='Malignant', ascending=False)

# Generate stacked bar plot
# 'tab20' colormap is used to distinguish multiple cell types
ax = df_sorted.plot(
    kind='bar', 
    stacked=True, 
    figsize=(max(15, len(df)//5), 8),  # Dynamically scale width based on sample count
    width=0.85, 
    colormap='tab20'
)

#-------------------------------------------------------------------------------
# 3. Plot Formatting and Customization
#-------------------------------------------------------------------------------

plt.title('Cell Composition Gradient (Sorted by Malignant Fraction)', fontsize=18)
plt.ylabel('Fraction (0.0 - 1.0)', fontsize=14)
plt.xlabel('Samples', fontsize=14)

# Place legend outside the plot for better readability
plt.legend(
    bbox_to_anchor=(1.02, 1), 
    loc='upper left', 
    title='Cell Types', 
    fontsize=10
)

# Adjust x-axis labels depending on number of samples
# Hide labels if too many samples to avoid clutter
if len(df_sorted) > 50:
    plt.xticks([])
    plt.xlabel(f'Samples (n={len(df_sorted)})', fontsize=14)
else:
    plt.xticks(rotation=90, fontsize=8)

plt.tight_layout()

#-------------------------------------------------------------------------------
# 4. Save the Figure
#-------------------------------------------------------------------------------

# Create output directory if it does not exist
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Save figure with high resolution
plt.savefig(os.path.join(save_dir, save_name), dpi=300)

print(f"Success! Plot saved to: {os.path.join(save_dir, save_name)}")
