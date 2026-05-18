"""
Script: Statistical Comparison of Statescope-Derived Malignant Fractions
        Between Group 3 and Group 4 Glioma Samples

Description:
This script evaluates differences in tumor purity, estimated as the malignant 
cell fraction derived from Statescope bulk RNA-seq deconvolution, between 
Group 3 and Group 4 glioma samples.

The script:
1. Loads Statescope deconvolution output and group classification data.
2. Merges datasets and filters for samples belonging to selected groups 
   (default: Group 3 and Group 4).
3. Extracts the malignant (tumor purity) fraction for each group.
4. Performs statistical testing to compare tumor purity between groups:
   - Assesses normality using the Shapiro–Wilk test.
   - Assesses homogeneity of variance using Levene’s test.
   - Automatically selects the appropriate statistical test:
       * Student’s t-test (equal variance)
       * Welch’s t-test (unequal variance)
       * Mann–Whitney U test (non-parametric)
   - Handles small sample sizes by avoiding unreliable normality assumptions.
5. Outputs summary statistics and detailed test results.
6. Visualizes group differences using a boxplot with overlaid individual sample points.
7. Saves the resulting figure for downstream analysis and reporting.

The resulting outputs provide insight into whether tumor purity differs 
between groups, which is important for interpreting downstream differences 
in tumor microenvironment composition.

Author: Mabel Pronk
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import shapiro, levene, mannwhitneyu, ttest_ind
from statannotations.Annotator import Annotator
#-------------------------------------------------------------------------------
# 1. Load and Filter Data
#-------------------------------------------------------------------------------
# Define file paths for deconvolution output and group classification
deconv_path = '/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output/fractions3.csv'
Group_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/final/classification_with_groups.csv'

# Load datasets (index = sample IDs)
df_deconv = pd.read_csv(deconv_path, index_col=0)
df_Groups = pd.read_csv(Group_path, index_col=0)

# Merge datasets and retain only samples in Group 3 and Group 4
df = df_deconv.join(df_Groups['Group']).dropna(subset=['Group'])
df_34 = df[df['Group'].isin(['Group 3', 'Group 4'])].copy()

#-------------------------------------------------------------------------------
# 2. Statistical Testing for 'Malignant' Fraction
#-------------------------------------------------------------------------------
# Extract malignant fractions for each group
g3_purity = df_34[df_34['Group'] == 'Group 3']['Malignant'].dropna()
g4_purity = df_34[df_34['Group'] == 'Group 4']['Malignant'].dropna()

n3, n4 = len(g3_purity), len(g4_purity)

# --- Normality Check (Shapiro-Wilk) ---
stat3, p_norm3 = shapiro(g3_purity) if n3 > 3 else (0, 0)
stat4, p_norm4 = shapiro(g4_purity) if n4 > 3 else (0, 0)
is_normal = (p_norm3 > 0.05 and p_norm4 > 0.05)

# --- Homogeneity of Variance (Levene's Test) ---
stat_l, p_lev = levene(g3_purity, g4_purity) if (n3 > 1 and n4 > 1) else (0, 0)
is_equal_var = (p_lev > 0.05)

# --- Select Appropriate Statistical Test ---
if is_normal and is_equal_var:
    test_name = "T-test"
    stat, p_val = ttest_ind(g3_purity, g4_purity, equal_var=True)
elif is_normal and not is_equal_var:
    test_name = "Welch's T"
    stat, p_val = ttest_ind(g3_purity, g4_purity, equal_var=False)
else:
    test_name = "Mann-Whitney"
    stat, p_val = mannwhitneyu(g3_purity, g4_purity, alternative='two-sided')

# Print summary of statistical results
print("--- Purity Comparison: Group 3 vs Group 4 ---")
print(f"Sample Sizes: Group 3 (n={n3}), Group 4 (n={n4})")
print(f"Normality (Shapiro p): G3={p_norm3:.4f}, G4={p_norm4:.4f}")
print(f"Variance (Levene p):  {p_lev:.4f}")
print(f"Statistical Test:    {test_name}")
print(f't-test: {stat}')
print(f"P-VALUE:             {p_val:.4f}")
print("-" * 45)

#-------------------------------------------------------------------------------
# 3. Plotting
#-------------------------------------------------------------------------------
# --- 3. Plotting (Refined for Consistency) ---
plt.figure(figsize=(7, 8))
sns.set_style("ticks")

ax = sns.boxplot(
    data=df_34, 
    x='Group', 
    y='Malignant', 
    order=['Group 3', 'Group 4'],
    palette={'Group 3': '#1f77b4', 'Group 4': '#ff7f0e'},
    width=0.5,
    linewidth=1.2
)

# Overlay individual sample points
sns.stripplot(data=df_34, x='Group', y='Malignant',order=['Group 3', 'Group 4'], color='black', alpha=0.3)

#pairs = [("Group 3", "Group 4")]
#annotator = Annotator(ax, pairs, data=df_34, x='Group', y='Malignant', order=['Group 3', 'Group 4'])
#annotator.configure(text_format="star", loc="inside", fontsize=12)
#annotator.set_pvalues([p_val]) # Uses the p_val calculated in Section 2
#annotator.annotate()

# Updated Title: T-stat and P-value in brackets
# Using 'stat' from the testing section above
plt.title(f'Tumor Purity \n({test_name}: t = {stat:.2f}, p = {p_val:.3f})', fontsize=15, pad=20)

plt.ylabel('Statescope malignant fraction estimate', fontsize=13)
plt.xlabel('Group', fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.ylim(0, 1.1)

# Black surrounding frame
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(1.2)

# Lighter Grid Lines
plt.grid(axis='y', linestyle='--', color='grey', alpha=0.3, linewidth=0.8)

plt.tight_layout()
#-------------------------------------------------------------------------------
# 4. Save Figure
#-------------------------------------------------------------------------------
save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output"
plt.savefig(os.path.join(save_dir, "purity_comparison_G3_G4.png"), dpi=300)
print(f"Purity boxplot saved to: {os.path.join(save_dir, 'purity_comparison_G3_G4.png')}")
