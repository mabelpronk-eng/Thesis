"""
Script: Comparison of Microenvironment Composition Between Groups

Description:
This script analyzes differences in tumor microenvironment (TME) composition 
between 2 groups of glioma samples using cell fractions derived from 
Statescope bulk RNA-seq deconvolution.

The script:
1. Loads deconvolution output and group classification data.
2. Merges datasets and filters for samples belonging to selected groups 
   (This can be adapted).
3. Excludes malignant and selected non-immune cell types to focus on the 
   microenvironment composition (this selection can be modified if needed).
4. Re-normalizes remaining cell-type fractions per sample to obtain relative 
   composition within the TME.
5. Visualizes differences in cell-type distributions using boxplots.
6. Performs statistical testing per cell type:
   - Assesses normality (Shapiro–Wilk test) and homogeneity of variance (Levene’s test).
   - Selects the appropriate test automatically:
       * Student’s t-test (equal variance)
       * Welch’s t-test (unequal variance)
       * Mann–Whitney U test (non-parametric)
   - Handles small sample sizes by skipping unreliable normality assessments.
7. Applies Benjamini–Hochberg (FDR) correction to account for multiple testing.
8. Outputs detailed statistical diagnostics and a summary table of raw and 
   adjusted p-values.
9. Saves the resulting visualization for downstream analysis.

The resulting figure and statistical summaries provide insight into differences 
in microenvironment composition between groups while controlling for multiple 
testing and statistical assumptions.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import shapiro, levene, mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests 
from statannotations.Annotator import Annotator

#-------------------------------------------------------------------------------
# 0. GROUP SELECTION (Select which groups to compare here)
#-------------------------------------------------------------------------------
# Change these two names to compare different groups (e.g., ['Group 1', 'Group 2'])
group_order = ['Group 1', 'Group 3']

g1_name = group_order[0]
g2_name = group_order[1]

# Define color palette logic
# Keeps Group 3/4 colors stable; assigns Group 1/2 different colors
master_colors = {
    'Group 1': '#2ca02c', # Green
    'Group 2': '#d62728', # Red
    'Group 3': '#1f77b4', # Blue
    'Group 4': '#ff7f0e'  # Orange
}
current_palette = {g: master_colors.get(g, '#7f7f7f') for g in group_order}

#-------------------------------------------------------------------------------
# 1. Load and Filter Data
#-------------------------------------------------------------------------------
# Define file paths for deconvolution output and group annotations
deconv_path = '/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output/statescope/fractions3.csv'
Group_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/IFNE_final_visual/classification_with_groups.csv'

# Load datasets (index = sample IDs)
df_deconv = pd.read_csv(deconv_path, index_col=0)
df_Groups = pd.read_csv(Group_path, index_col=0)

df = df_deconv.join(df_Groups['Group']).dropna(subset=['Group'])

# Subset to the selected groups
df_subset = df[df['Group'].isin(group_order)].copy()

# Get sample counts dynamically
n_g1 = len(df_subset[df_subset['Group'] == g1_name])
n_g2 = len(df_subset[df_subset['Group'] == g2_name])

#-------------------------------------------------------------------------------
# 2. Exclude Malignant Cells and Re-normalize
#-------------------------------------------------------------------------------
# Remove cell that are not immune cells
df_micro = df_subset.drop(columns=['Malignant', 'Oligodendrocyte', 'Endothelial', 'Pericyte', 'Fibroblast'])

# Identify all remaining cell-type columns (exclude group label)
cell_types = [c for c in df_micro.columns if c != 'Group']

# Re-normalize fractions so that remaining cell types sum to 1 per sample
# This converts values into relative composition within the microenvironment
df_micro[cell_types] = df_micro[cell_types].div(df_micro[cell_types].sum(axis=1), axis=0)

# Convert to long format for plotting with seaborn
df_melted = df_micro.reset_index().melt(
    id_vars=['index', 'Group'], 
    value_vars=cell_types,
    var_name='Cell Type', 
    value_name='Relative Fraction'
)

#-------------------------------------------------------------------------------
# 3. Statistical Testing & BH Correction
#-------------------------------------------------------------------------------

raw_p_values = []
test_list = []
insufficient_data = [] # Track small sample sizes

print(f"\n{'Cell Type':<18} | {'Shapiro (p)':<12} | {'Levene (p)':<10} | {'Test Used':<15} | {'p-val (Raw)'}")
print("-" * 85)

for cell in cell_types:
    data1 = df_micro[df_micro['Group'] == g1_name][cell].dropna()
    data2 = df_micro[df_micro['Group'] == g2_name][cell].dropna()
    
    n1, n2 = len(data1), len(data2)

    if n1 <= 3 or n2 <= 3:
        insufficient_data.append(f"{cell} ({g1_name}: n={n1}, {g2_name}: n={n2})")
        p_norm1, p_norm2 = 0, 0 
        is_normal = False
    else:
        p_norm1 = shapiro(data1)[1]
        p_norm2 = shapiro(data2)[1]
        is_normal = (p_norm1 > 0.05 and p_norm2 > 0.05)

    stat_l, p_lev = levene(data1, data2) if (n1 > 1 and n2 > 1) else (0, 0)
    is_equal_var = (p_lev > 0.05)

    if is_normal and is_equal_var:
        test_name = "T-test"
        _, p_val = ttest_ind(data1, data2, equal_var=True)
    elif is_normal and not is_equal_var:
        test_name = "Welch's T"
        _, p_val = ttest_ind(data1, data2, equal_var=False)
    else:
        test_name = "Mann-Whitney"
        _, p_val = mannwhitneyu(data1, data2, alternative='two-sided')

    raw_p_values.append(p_val)
    test_list.append(test_name)
    shapiro_str = f"{p_norm1:.3f}/{p_norm2:.3f}"
    print(f"{cell:<18} | {shapiro_str:<12} | {p_lev:<10.3f} | {test_name:<15} | {p_val:.4f}")

_, adj_p_values, _, _ = multipletests(raw_p_values, method='fdr_bh')

#Overview table "Before vs After"
# --- 6. Final Summaries (RESTORED) ---
if insufficient_data:
    print(f"\nATTENTION: Normality testing skipped for {len(insufficient_data)} cell types (n <= 3):")
    for item in insufficient_data:
        print(f" - {item}")

print(f"\n{'='*85}")
print(f"{'CELL TYPE OVERVIEW: RAW P vs ADJUSTED FDR (BH)':^85}")
print(f"{'='*85}")
print(f"{'Cell Type':<18} | {f'n ({g1_name}/{g2_name})':<12} | {'Raw P':<10} | {'Adj Q (FDR)':<12} | {'Sig?'}")
print(f"{'-'*85}")

for i, cell in enumerate(cell_types):
    sig = "YES! *" if adj_p_values[i] < 0.05 else "no"
    n_info = f"{n_g1}/{n_g2}"
    print(f"{cell:<18} | {n_info:<12} | {raw_p_values[i]:<10.4f} | {adj_p_values[i]:<12.4f} | {sig}")
print(f"{'='*85}")

#-------------------------------------------------------------------------------
# 4. Plotting
#-------------------------------------------------------------------------------
plt.figure(figsize=(16, 8))
sns.set_style("ticks")

ax = sns.boxplot(
    data=df_melted, 
    x='Cell Type', 
    y='Relative Fraction', 
    hue='Group',
    hue_order=group_order,  # NEW: Forces Group 3 to the left
    palette=current_palette,
    linewidth=1.2,
    showfliers=False
)

# --- ADDED ANNOTATOR BLOCK ---
pairs = [((cell, g1_name), (cell, g2_name)) for cell in cell_types]
annotator = Annotator(ax, pairs, data=df_melted, x='Cell Type', y='Relative Fraction', hue='Group', hue_order=group_order)
annotator.configure(text_format="star", loc="inside", fontsize=12)
annotator.set_pvalues(adj_p_values) # Uses the FDR-adjusted values
annotator.annotate()

plt.title(f'Microenvironment Composition: {g1_name} vs {g2_name}', fontsize=16)
plt.ylabel('Relative Fraction')
plt.xlabel('Cell Type')
plt.xticks(rotation=45, ha='right')
plt.legend(title='Group', loc='upper right')


plt.tight_layout()

#-------------------------------------------------------------------------------
# 4. Save Figure
#-------------------------------------------------------------------------------
# Define output directory and filename
save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output/IFNE_CDKN2AB/"
save_name = f"boxplot_TME_{g1_name.replace(' ', '')}_vs_{g2_name.replace(' ', '')}.png"
save_path = os.path.join(save_dir, save_name)

if not os.path.exists(save_dir): os.makedirs(save_dir)
plt.savefig(save_path, dpi=300)
print(f"\nPlot saved to: {save_path}")



