"""
Script: Comparison of Microenvironment Composition Between Group 3 and Group 4

Description:
This script analyzes differences in tumor microenvironment (TME) composition 
between Group 3 and Group 4 glioma samples using cell fractions derived from 
Statescope bulk RNA-seq deconvolution.

The script:
1. Loads deconvolution output and group classification data.
2. Merges datasets and filters for samples belonging to selected groups 
   (default: Group 3 and Group 4, but this can be adapted).
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

#-------------------------------------------------------------------------------
# 1. Load and Filter Data
#-------------------------------------------------------------------------------
# Define file paths for deconvolution output and group annotations
deconv_path = '/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output/fractions3.csv'
Group_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/final/classification_with_groups.csv'

# Load datasets (index = sample IDs)
df_deconv = pd.read_csv(deconv_path, index_col=0)
df_Groups = pd.read_csv(Group_path, index_col=0)

# Merge datasets and retain only samples with group labels
df = df_deconv.join(df_Groups['Group']).dropna(subset=['Group'])

# Explicitly define order
# Subset to Groups 3 and 4 for comparison
group_order = ['Group 3', 'Group 4']
df_34 = df[df['Group'].isin(group_order)].copy()

# Get sample counts for the legend
n_g3 = len(df_34[df_34['Group'] == 'Group 3'])
n_g4 = len(df_34[df_34['Group'] == 'Group 4'])
#-------------------------------------------------------------------------------
# 2. Exclude Malignant Cells and Re-normalize
#-------------------------------------------------------------------------------
# Remove cell that are not immune cells
df_micro = df_34.drop(columns=['Malignant', 'Oligodendrocyte', 'Endothelial', 'Pericyte', 'Fibroblast'])

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
# 3. Plotting
#-------------------------------------------------------------------------------
plt.figure(figsize=(14, 8))
sns.set_style("ticks")

ax = sns.boxplot(
    data=df_melted, 
    x='Cell Type', 
    y='Relative Fraction', 
    hue='Group',
    hue_order=group_order,  # NEW: Forces Group 3 to the left
    palette={'Group 3': '#1f77b4', 'Group 4': '#ff7f0e'},
    linewidth=1.2,
    showfliers=True
)

plt.title('Microenvironment Composition: Group 3 vs Group 4', fontsize=15, pad=20)
plt.ylabel('Relative Fraction', fontsize=13)
plt.xlabel('Cell Type', fontsize=13)
plt.xticks(rotation=45, ha='right', fontsize=11)

# Black surrounding frame (Spines)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(1.2)

# Lighter Grid Lines (Back to the original subtle look)
plt.grid(axis='y', linestyle='--', color='grey', alpha=0.3, linewidth=0.8)
plt.legend(title="Group", loc='upper right', frameon=True, edgecolor='black', fontsize=11)
plt.tight_layout()

#-------------------------------------------------------------------------------
# 4. Save Figure
#-------------------------------------------------------------------------------
# Define output directory and filename
save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output"
save_path = os.path.join(save_dir, "boxplot_microenvironment_G3_G4_renormalized.png")
if not os.path.exists(save_dir): os.makedirs(save_dir)
plt.savefig(save_path, dpi=300)

#-------------------------------------------------------------------------------
# 5. Statistical Testing & BH Correction
#-------------------------------------------------------------------------------

raw_p_values = []
test_list = []
insufficient_data = [] # Track small sample sizes

print(f"\n{'Cell Type':<18} | {'Shapiro (p)':<12} | {'Levene (p)':<10} | {'Test Used':<15} | {'p-val (Raw)'}")
print("-" * 85)

for cell in cell_types:
    g3_data = df_micro[df_micro['Group'] == 'Group 3'][cell].dropna()
    g4_data = df_micro[df_micro['Group'] == 'Group 4'][cell].dropna()
    
     # --- A. Normality Test (Shapiro-Wilk) ---
     # Check if both groups follow a normal distribution
    n3, n4 = len(g3_data), len(g4_data)

    if n3 <= 3 or n4 <= 3:
        insufficient_data.append(f"{cell} (G3: n={n3}, G4: n={n4})")
        p_norm3, p_norm4 = 0, 0 
        is_normal = False
    else:
        p_norm3 = shapiro(g3_data)[1]
        p_norm4 = shapiro(g4_data)[1]
        is_normal = (p_norm3 > 0.05 and p_norm4 > 0.05)

    # --- B. Homogeneity of Variance (Levene’s test) ---
    stat_l, p_lev = levene(g3_data, g4_data) if (n3 > 1 and n4 > 1) else (0, 0)
    is_equal_var = (p_lev > 0.05)

    if is_normal and is_equal_var:
        test_name = "T-test"
        _, p_val = ttest_ind(g3_data, g4_data, equal_var=True)
    elif is_normal and not is_equal_var:
        test_name = "Welch's T"
        _, p_val = ttest_ind(g3_data, g4_data, equal_var=False)
    else:
        test_name = "Mann-Whitney"
        _, p_val = mannwhitneyu(g3_data, g4_data, alternative='two-sided')

    raw_p_values.append(p_val)
    test_list.append(test_name)
    shapiro_str = f"{p_norm3:.3f}/{p_norm4:.3f}"
    print(f"{cell:<18} | {shapiro_str:<12} | {p_lev:<10.3f} | {test_name:<15} | {p_val:.4f}")

#Apply Benjamini-Hochberg Correction
_, adj_p_values, _, _ = multipletests(raw_p_values, method='fdr_bh')

#Overview table "Before vs After"
# --- 6. Final Summaries ---
if insufficient_data:
    print(f"\nATTENTION: Normality testing skipped for {len(insufficient_data)} cell types (n <= 3):")
    for item in insufficient_data:
        print(f" - {item}")

print(f"\n{'='*85}")
print(f"{'CELL TYPE OVERVIEW: RAW P vs ADJUSTED FDR (BH)':^85}")
print(f"{'='*85}")
print(f"{'Cell Type':<18} | {'n (G3/G4)':<12} | {'Raw P':<10} | {'Adj Q (FDR)':<12} | {'Sig?'}")
print(f"{'-'*85}")

for i, cell in enumerate(cell_types):
    sig = "YES! *" if adj_p_values[i] < 0.05 else "no"
    n_info = f"{n_g3}/{n_g4}"
    print(f"{cell:<18} | {n_info:<12} | {raw_p_values[i]:<10.4f} | {adj_p_values[i]:<12.4f} | {sig}")
print(f"{'='*85}")
