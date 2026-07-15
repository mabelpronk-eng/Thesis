"""
============================================================
Global Mean Cell-Type Composition (Stacked Barplot)
============================================================

Description:
This script computes and visualizes the average cell-type
composition across a TCGA cohort based on deconvolution
outputs (Statescope or CIBERSORTx).

Steps:
1. Load deconvolution fraction matrix
2. Compute mean cell-type fractions across all samples
3. Sort cell types by abundance
4. Map predefined cell-type colors (from JSON file)
5. Generate a stacked barplot of global mean composition
6. Save the figure for downstream analysis

Output:
- Stacked barplot showing mean cohort composition
- Consistent cell-type coloring across all figures

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import json

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,

    "font.family": "Times New Roman",
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ==================================================
# Load predefined cell-type color mapping
# ==================================================
color_path = '/net/beegfs/users/P086608/celltype_colors.json'
with open(color_path, "r") as f:
    CELLTYPE_COLORS = json.load(f)

# ==================================================
# 1. LOAD DECONVOLUTION DATA
# ==================================================
# Choose between Statescope or CIBERSORTx output
deconv_path = '/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output/statescope/fractions3.csv'
#deconv_path = '/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output/cibersort_fractions_TCGAbulk_gbm.tsv'

# Output directory for figure
#save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output"
save_dir = '/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output'
save_name = "global_mean_composition_stacked.pdf"

# Read deconvolution fraction matrix (samples × cell types)
df_deconv = pd.read_csv(deconv_path,  index_col=0)

# ==================================================
# 2. COMPUTE GLOBAL MEAN COMPOSITION
# ==================================================
# Calculate the mean of every column (cell type) across all samples
# Then convert the resulting Series into a single-row DataFrame
mean_series = df_deconv.mean()

# Sort by abundance so the legend and stack order are logical (largest at bottom)
mean_series = mean_series.sort_values(ascending=False)

# Convert to single-row DataFrame for plotting
df_global_mean = mean_series.to_frame().T
df_global_mean.index = ['All Samples (N={})'.format(len(df_deconv))]

# ==================================================
# 3. MAP CELL TYPES TO CONSISTENT COLORS
# ==================================================
cell_types = df_global_mean.columns
# Assign predefined colors; fallback to grey if missing
colors = [CELLTYPE_COLORS.get(ct, "#808080") for ct in cell_types]

# ==================================================
# 4. GENERATE STACKED BARPLOT
# ==================================================
fig, ax = plt.subplots(figsize=(6, 7))

df_global_mean.plot(
    kind='bar', 
    stacked=True, 
    ax=ax,
    color=colors, 
    width=0.5,
    edgecolor='white',
    linewidth=0.5
)

# ==================================================
# 5. FORMATTING
# ==================================================
plt.title('Stacked barplot of average cohort composition', fontsize=14, pad=20)
plt.ylabel('Mean Fraction', fontsize=14)
plt.xlabel('') # No x-axis label (single bar representation)
plt.xticks(rotation=0, fontsize=12) 
plt.ylim(0, 1.0) 

# Legend: Reversed to match the stacking order (Top of bar = Top of legend)
handles, labels = ax.get_legend_handles_labels()
plt.legend(
    reversed(handles), reversed(labels),
    title='Cell Types', 
    bbox_to_anchor=(1.05, 1), 
    loc='upper left', 
    fontsize=12,
    frameon=False
)

# ==================================================
# 6. SAVE FIGURE
# ==================================================
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.tight_layout()
plt.savefig(os.path.join(save_dir, save_name), bbox_inches='tight')
print(f"Success! Global mean stacked plot saved to: {os.path.join(save_dir, save_name)}")
plt.show()
