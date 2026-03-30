"""
Script: Comparison of Microenvironment Composition Between Group 3 and Group 4

Description:
This script analyzes differences in tumor microenvironment composition between 
Group 3 and Group 4 glioma samples using cell fractions derived from Statescope 
bulk RNA-seq deconvolution.

The script:
1. Loads deconvolution output and group classification data.
2. Merges datasets and filters for samples belonging to Group 3 and Group 4. You can also look at other groups if preferred. 
3. Excludes malignant cell fractions to focus on the microenvironment. You can also exclude other cells if preferred.
4. Re-normalizes remaining cell-type fractions per sample to obtain relative composition.
5. Visualizes differences in cell-type distributions using boxplots.
6. Performs statistical testing per cell type:
   - Assesses normality (Shapiro–Wilk test) and homogeneity of variance (Levene’s test).
   - Applies either a t-test (standard or Welch’s) or Mann–Whitney U test as appropriate.
7. Outputs statistical results and saves the visualization for downstream analysis.

The resulting figure and statistical summaries provide insight into differences in 
microenvironment composition between the two groups.

Author: Mabel Pronk
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu

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

# Subset to Groups 3 and 4 for comparison
df_34 = df[df['Group'].isin(['Group 3', 'Group 4'])].copy()

#-------------------------------------------------------------------------------
# 2. Exclude Malignant Cells and Re-normalize
#-------------------------------------------------------------------------------
# Remove malignant cell fraction to focus on the microenvironment
df_micro = df_34.drop(columns=['Malignant'])

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
# Initialize figure
plt.figure(figsize=(14, 7))

# Create boxplot comparing cell-type distributions between groups
ax = sns.boxplot(
    data=df_melted, 
    x='Cell Type', 
    y='Relative Fraction', 
    hue='Group',
    palette={'Group 3': '#1f77b4', 'Group 4': '#ff7f0e'},  # consistent color scheme
    showfliers=True  # retain outliers
)

# Customize plot appearance
plt.title('Microenvironment Composition: Group 3 vs Group 4\n(Normalized: Malignant Cells Excluded)', fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Relative Fraction of Microenvironment')
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()

#-------------------------------------------------------------------------------
# 4. Save Figure
#-------------------------------------------------------------------------------
# Define output directory and filename
save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output"
save_path = os.path.join(save_dir, "boxplot_microenvironment_G3_G4_renormalized.png")

# Create directory if it does not exist
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Save figure to file
plt.savefig(save_path, dpi=300)
print(f"Boxplot saved to: {save_path}")

#-------------------------------------------------------------------------------
# 5. Statistical Testing
#-------------------------------------------------------------------------------
from scipy.stats import shapiro, levene, mannwhitneyu, ttest_ind

# Print number of samples per group
group_counts = df_micro['Group'].value_counts()
print("--- Group Sample Sizes ---")
for grp, count in group_counts.items():
    print(f"{grp}: n = {count}")
print("-" * 30)

# Prepare storage for results
results = []
insufficient_data = []  # track cell types with very small sample sizes

# Identify cell-type columns again
cell_types = [c for c in df_micro.columns if c != 'Group']

# Print header for results table
print(f"{'Cell Type':<18} | {'Shapiro (p)':<12} | {'Levene (p)':<10} | {'Test Used':<15} | {'p-val (Raw)'}")
print("-" * 85)

# Loop through each cell type and perform statistical testing
for cell in cell_types:
    # Extract data per group
    g3_data = df_micro[df_micro['Group'] == 'Group 3'][cell].dropna()
    g4_data = df_micro[df_micro['Group'] == 'Group 4'][cell].dropna()
    
    # --- A. Normality Test (Shapiro-Wilk) ---
    # Check if both groups follow a normal distribution
    n3, n4 = len(g3_data), len(g4_data)

    # Skip normality testing if sample size is too small (n <= 3)
    if n3 <= 3 or n4 <= 3:
        insufficient_data.append(f"{cell} (G3: n={n3}, G4: n={n4})")
        p_norm3, p_norm4 = 0, 0 
        is_normal = False
    else:
        stat3, p_norm3 = shapiro(g3_data)
        stat4, p_norm4 = shapiro(g4_data)
        is_normal = (p_norm3 > 0.05 and p_norm4 > 0.05)
    
    # --- B. Homogeneity of Variance (Levene’s test) ---
    stat_l, p_lev = levene(g3_data, g4_data) if (len(g3_data) > 1 and len(g4_data) > 1) else (0, 0)
    is_equal_var = (p_lev > 0.05)

    # --- C. Select Appropriate Statistical Test ---
    if is_normal and is_equal_var:
        # Standard independent t-test
        test_name = "T-test"
        stat, p_val = ttest_ind(g3_data, g4_data, equal_var=True)
    elif is_normal and not is_equal_var:
        # Welch’s t-test (unequal variances)
        test_name = "Welch's T"
        stat, p_val = ttest_ind(g3_data, g4_data, equal_var=False)
    else:
        # Mann–Whitney U test (non-parametric)
        test_name = "Mann-Whitney"
        stat, p_val = mannwhitneyu(g3_data, g4_data, alternative='two-sided')

    # Print summary for each cell type
    shapiro_str = f"{p_norm3:.3f}/{p_norm4:.3f}"
    print(f"{cell:<18} | {shapiro_str:<12} | {p_lev:<10.3f} | {test_name:<15} | {p_val:.4f}")
    
    # Store results
    results.append({
        'Cell Type': cell,
        'p_raw': p_val,
        'Test': test_name,
        'n_G3': len(g3_data),
        'n_G4': len(g4_data)
    })

# Report cell types with insufficient sample size
if insufficient_data:
    print(f"ATTENTION: Normality testing skipped for {len(insufficient_data)} cell types (n <= 3):")
    for item in insufficient_data:
        print(f" - {item}")
else:
    print("All cell types had sufficient samples (n > 3) for normality testing.")
