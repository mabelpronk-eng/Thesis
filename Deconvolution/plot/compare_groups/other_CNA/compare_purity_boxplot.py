"""
============================================================
Script: Tumor Purity Comparison Across CNA Status Groups
============================================================

Description:
This script evaluates differences in tumor purity (Malignant fraction)
between samples with different copy number alteration (CNA) statuses
(e.g., PTEN deletion vs neutral).

Workflow:
1. Load GISTIC-derived CNA status annotations.
2. Load Statescope deconvolution results (cell-type fractions).
3. Harmonize sample IDs between datasets.
4. Merge malignant fraction with CNA status information.
5. Filter samples for selected comparison groups.
6. Perform statistical testing on tumor purity:
   - Shapiro–Wilk test for normality
   - Levene’s test for variance equality
   - Automatic selection of:
       * Student’s t-test
       * Welch’s t-test
       * Mann–Whitney U test
7. Visualize differences using boxplots with:
   - Individual sample overlay (stripplot)
   - Significance annotation (statannotations)
8. Save publication-quality figure.

Output:
- Boxplot comparing Malignant fraction across CNA groups
- Statistical test summary printed to console
- Saved figure in PDF format

============================================================
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import shapiro, levene, ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator

# ============================================================
# Global plotting configuration (publication-ready styling)
# ============================================================
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


# ============================================================
# 1. INPUT PATHS AND PARAMETERS
# ============================================================
gistic_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/GISTIC/other_CNA/Focal_Aberrations_with_Status.xlsx'
deconv_path = '/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output/statescope/fractions3.csv'
save_dir = "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output/Other_focal/PTEN/"

target_gene_status = "PTEN_Status" 

# Define groups to compare here
comparison_groups = ['Neutral', 'Deletion']

# ============================================================
# 2. LOAD AND MERGE DATASETS
# ============================================================
df_gistic = pd.read_excel(gistic_path)
df_deconv = pd.read_csv(deconv_path, index_col=0)

# Harmonize sample IDs to first 12 characters (TCGA format)
df_gistic['Match_ID'] = df_gistic['Sample'].str.slice(0, 12)
df_deconv.index = df_deconv.index.str.slice(0, 12)

# Merge CNA status with malignant fraction
df = df_deconv[['Malignant']].join(df_gistic.set_index('Match_ID')[[target_gene_status]], how='inner')

# ============================================================
# 3. FILTER SELECTED GROUPS
# ============================================================
df_2group = df[df[target_gene_status].isin(comparison_groups)].copy()

# Output sample counts
print(f"\n--- Sample Counts for {target_gene_status} ---")
print(df_2group[target_gene_status].value_counts())
print("-" * 40)

# ============================================================
# 4. STATISTICAL TESTING (MALIGNANT FRACTION)
# ============================================================
# Dynamically assign data based on comparison_groups list
g1_name = comparison_groups[0]
g2_name = comparison_groups[1]

g1_purity = df_2group[df_2group[target_gene_status] == g1_name]['Malignant'].dropna()
g2_purity = df_2group[df_2group[target_gene_status] == g2_name]['Malignant'].dropna()

n1, n2 = len(g1_purity), len(g2_purity)

# --- Normality Check (Shapiro-Wilk) ---
_, p_norm1 = shapiro(g1_purity) if n1 >= 3 else (0, 0)
_, p_norm2 = shapiro(g2_purity) if n2 >= 3 else (0, 0)
is_normal = (p_norm1 > 0.05 and p_norm2 > 0.05)

# --- Homogeneity of Variance (Levene's Test) ---
_, p_lev = levene(g1_purity, g2_purity) if (n1 > 1 and n2 > 1) else (0, 0)
is_equal_var = (p_lev > 0.05)

# --- Select Appropriate Statistical Test ---
if is_normal and is_equal_var:
    test_name = "T-test"
    stat, p_val = ttest_ind(g1_purity, g2_purity, equal_var=True)
elif is_normal and not is_equal_var:
    test_name = "Welch's T"
    stat, p_val = ttest_ind(g1_purity, g2_purity, equal_var=False)
else:
    test_name = "Mann-Whitney"
    stat, p_val = mannwhitneyu(g1_purity, g2_purity, alternative='two-sided')

# Print summary
print(f"--- Purity Comparison: {g1_name} vs {g2_name} ({target_gene_status}) ---")
print(f"Sample Sizes: {g1_name} (n={n1}), {g2_name} (n={n2})")
print(f"Normality (Shapiro p): {g1_name}={p_norm1:.4f}, {g2_name}={p_norm2:.4f}")
print(f"Variance (Levene p):  {p_lev:.4f}")
print(f"Statistical Test Used: {test_name}")
print(f"P-VALUE: {p_val:.4f}")
print(f"STATISTIC ({test_name}): {stat:.4f}")
print("-" * 60)

# ============================================================
# 5. VISUALIZATION
# ============================================================
plt.figure(figsize=(4, 6))

# Define a palette mapping to ensure colors stay consistent
palette_map = {'Neutral': '#bdc3c7', 'Amplification': '#e74c3c', 'Deletion': '#3498db'}

ax = sns.boxplot(
    data=df_2group, 
    x=target_gene_status, 
    y='Malignant', 
    order=comparison_groups,
    palette=palette_map,
    width=0.5,
    linewidth=1.2
)

# Overlay individual sample points
sns.stripplot(data=df_2group, x=target_gene_status, y='Malignant', order=comparison_groups, color='black', alpha=0.3)

# Add Significance Line
pairs = [(comparison_groups[0], comparison_groups[1])]
annotator = Annotator(ax=plt.gca(), pairs=pairs, data=df_2group, 
                      x=target_gene_status, y='Malignant', order=comparison_groups)

# Configure the look (uses stars: * for p<0.05, ** for p<0.01, etc.)
annotator.configure(text_format="star", loc="inside")
annotator.set_pvalues([p_val])
annotator.annotate()

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(1.1)

# Customize plot appearance
plt.title(f'Tumor Purity Distribution ({target_gene_status}) \n({test_name}: t = {stat:.2f}, p = {p_val:.3f})', pad=20)
plt.ylabel('Statescope Malignant Fraction')
plt.ylim(0, 1.1) 
plt.grid(axis='y', linestyle='--', alpha=0.3, linewidth=0.8)
plt.tight_layout()

# ============================================================
# 6. SAVE OUTPUT
# ============================================================
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(os.path.join(save_dir, f"Purity_Stats_{g1_name}_vs_{g2_name}_{target_gene_status}.pdf"))
print(f"Purity boxplot saved to: {os.path.join(save_dir, f'Purity_Stats_{g1_name}_vs_{g2_name}_{target_gene_status}.png')}")
