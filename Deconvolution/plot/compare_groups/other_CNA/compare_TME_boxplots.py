# ============================================================
# Script: TME Composition Analysis by Copy Number Status 
# ============================================================
#
# Description:
# This script integrates bulk tumor deconvolution data with GISTIC-derived
# copy number status to investigate how genomic alterations (e.g., CDK4 amplification)
# are associated with changes in tumor microenvironment (TME) composition.
#
# Workflow:
# 1. Load:
#    - Gene-level copy number status (GISTIC, with categorical annotations)
#    - Cell type fractions (deconvolution output, e.g. Statescope)
# 2. Harmonize sample IDs and merge datasets.
# 3. Optionally exclude specific cell types (e.g., malignant or non-immune)
#    and re-normalize remaining fractions.
# 4. Subset samples into comparison groups (e.g., Neutral vs Amplification).
# 5. Perform statistical testing per cell type:
#    - Normality (Shapiro-Wilk)
#    - Variance equality (Levene’s test)
#    - Group differences (Mann–Whitney U test)
# 6. Apply multiple testing correction (Benjamini–Hochberg FDR).
# 7. Visualize results using boxplots with statistical annotations.
# 8. Print a summary table comparing raw and adjusted p-values.
#
# Outputs:
# - Boxplot showing differences in cell type fractions between CNA groups
# - Console output with detailed statistical diagnostics
# - Summary table indicating significance before and after correction
#
# Purpose:
# To identify associations between focal genomic aberrations (e.g., CDK4 amplification)
# and shifts in tumor microenvironment composition, supporting biological interpretation
# of CNA-driven TME interactions.
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import shapiro, levene, mannwhitneyu
from statannotations.Annotator import Annotator

# ------------------------------------------------------------
# Plotting configuration (publication-ready styling)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 1. INPUT FILES AND PARAMETERS
# ------------------------------------------------------------
# GISTIC-derived CNA annotations (sample-level gene status)
gistic_with_groups_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/GISTIC/other_CNA/Focal_Aberrations_with_Status.xlsx'

# Deconvolution output (cell type fractions per sample)
deconv_path = '/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output/statescope/fractions3.csv'
#deconv_path = '/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output/cibersort_fractions_TCGAbulk_gbm.tsv'

# Output directory for plots
save_dir = "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output/Other_focal/PTEN"
#save_dir = '/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output'

# CNA status column of interest
target_gene_status = "PTEN_Status" 

# Groups to compare (must match values in GISTIC table)
comparison_groups = ['Neutral', 'Deletion'] 
g1_name, g2_name = comparison_groups[0], comparison_groups[1]

# ------------------------------------------------------------
# 2. ANALYSIS CONFIGURATION (CELL TYPE FILTERING)
# ------------------------------------------------------------
# Defines which cell types are excluded from analysis
# depending on biological focus
# Options: 'all' (no exclusion), 'no_malignant' (only Malignant out), 'immune_only' (all non-immune out)
analysis_mode = 'immune_only' 

# Define exclusion rules based on the mode
exclusion_rules = {
    'all': [],
    'no_malignant': ['Malignant'],
    'immune_only': ['Malignant', 'Oligodendrocyte', 'Endothelial', 'Pericyte', 'Fibroblast']
}
to_exclude = exclusion_rules[analysis_mode]

# ------------------------------------------------------------
# 3. DATA LOADING AND SAMPLE ALIGNMENT
# ------------------------------------------------------------
# Load CNA annotations and deconvolution results
df_gistic = pd.read_excel(gistic_with_groups_path)
df_deconv = pd.read_csv(deconv_path, index_col=0)

df_gistic['Match_ID'] = df_gistic['Sample'].str.slice(0, 12)
df_deconv.index = df_deconv.index.str.slice(0, 12)

df = df_deconv.join(df_gistic.set_index('Match_ID')[[target_gene_status]], how='inner')

# ------------------------------------------------------------
# 4. DATA QUALITY CHECK (GROUP COUNTS)
# ------------------------------------------------------------
df_filtered_audit = df[df[target_gene_status].isin(comparison_groups)]
print(f"\n--- Analysis Mode: {analysis_mode} (Excluding: {to_exclude}) ---")
print(df_filtered_audit[target_gene_status].value_counts())
print("-" * 45)

# ------------------------------------------------------------
# 5. PREPROCESSING (EXCLUSION + RENORMALIZATION)
# ------------------------------------------------------------
# Remove selected cell types depending on analysis mode
# We only drop columns that actually exist in the dataframe to prevent errors
df_micro = df.drop(columns=[col for col in to_exclude if col in df.columns])

# Identify remaining cell types
cell_types = [c for c in df_micro.columns if c != target_gene_status]

# Renormalize cell fractions after exclusion
# (ensures compositional structure is preserved)
df_micro[cell_types] = df_micro[cell_types].div(df_micro[cell_types].sum(axis=1), axis=0)

# Prepare long-format dataframe for plotting
df_plot = df_micro[df_micro[target_gene_status].isin(comparison_groups)]
df_melted = df_plot.reset_index().melt(id_vars=['index', target_gene_status], value_vars=cell_types)

# ------------------------------------------------------------
# 6. STATISTICAL TESTING PER CELL TYPE
# ----------------------------------------------------------
p_values_dict = {}

# Color mapping for CNA groups (plot consistency)
palette_map = {'Neutral': '#bdc3c7', 'Amplification': '#e74c3c', 'Deletion': '#3498db'}

print(f"\n{'='*105}")
print(f" STATISTICAL DIAGNOSTICS: {g1_name} vs {g2_name}")
print(f"{'='*105}")
print(f"{'Cell Type':<18} | {'N (G1/G2)':<12} | {f'Shap {g1_name[:3]}':<12} | {f'Shap {g2_name[:3]}':<12} | {'Levene p':<10} | {'MWU p'}")
print(f"{'-'*105}")

# Loop over cell types and perform statistical tests
for cell in cell_types:
    g1_data = df_micro[df_micro[target_gene_status] == g1_name][cell].dropna()
    g2_data = df_micro[df_micro[target_gene_status] == g2_name][cell].dropna()
    
    n1, n2 = len(g1_data), len(g2_data)

    # Only test if sample size is sufficient for normality assumptions
    if n1 >= 3 and n2 >= 3:
        _, p_shap1 = shapiro(g1_data)
        _, p_shap2 = shapiro(g2_data)
        _, p_levene = levene(g1_data, g2_data)
        _, p_mwu = mannwhitneyu(g1_data, g2_data, alternative='two-sided')
        p_values_dict[cell] = p_mwu
        print(f"{cell:<18} | {n1:<3}/{n2:<6} | {p_shap1:<12.4f} | {p_shap2:<12.4f} | {p_levene:<10.4f} | {p_mwu:.4f}")
    else:
        # fallback for low sample sizes
        p_values_dict[cell] = 1.0
        print(f"{cell:<18} | {n1:<3}/{n2:<6} | {'Low N':<12} | {'Low N':<12} | {'N/A':<10} | 1.0000")

# ------------------------------------------------------------
# 7. MULTIPLE TESTING CORRECTION (FDR)
# -----------------------------------------------------------
from statsmodels.stats.multitest import multipletests

# Extract the raw p-values in the order of cell_types
raw_p_values = [p_values_dict[cell] for cell in cell_types]

# Apply Benjamini-Hochberg (fdr_bh) correction
# This returns: [is_rejected, corrected_p_values, alphacSidac, alphacBonf]
_, corrected_p_values, _, _ = multipletests(raw_p_values, method='fdr_bh')

# Update the dictionary with corrected values to use in the plot
corrected_p_dict = dict(zip(cell_types, corrected_p_values))

# ------------------------------------------------------------
# 8. VISUALIZATION (BOXPLOT + SIGNIFICANCE ANNOTATION)
# ----------------------------------------------------------
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=df_melted, 
    x='variable', 
    y='value', 
    hue=target_gene_status, 
    hue_order=comparison_groups,
    palette=palette_map, 
    linewidth=1.2,
    showfliers=True
)

ax.yaxis.grid(True, linestyle='--', linewidth=0.7, alpha=0.6)

pairs = [((cell, g1_name), (cell, g2_name)) for cell in cell_types]
annotator = Annotator(ax, pairs, data=df_melted, x='variable', y='value', hue=target_gene_status, hue_order=comparison_groups)
annotator.configure(text_format="star", loc="inside")
annotator.set_pvalues([corrected_p_dict[cell] for cell in cell_types])
annotator.annotate()

# Black surrounding frame
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(1.1)

plt.title(f'Microenvironment Composition ({target_gene_status}): {g1_name} vs {g2_name}')
plt.xticks(rotation=45, ha='right')
plt.xlabel('Cell Type')
plt.ylabel('Relative Fraction')
plt.legend(title='Group', loc='upper right')
plt.tight_layout()

# ------------------------------------------------------------
# 9. SAVE OUTPUT
# -----------------------------------------------------------
if not os.path.exists(save_dir): os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f"TME_{analysis_mode}_{g1_name}_vs_{g2_name}_{target_gene_status}.pdf"))

# ------------------------------------------------------------
# 10. SUMMARY TABLE (RAW VS CORRECTED P-VALUES)
# ----------------------------------------------------------
print(f"\n{'='*95}")
print(f"{'CELL TYPE':<20} | {'RAW P-VAL':<12} | {'ADJ Q-VAL':<12} | {'STATUS'}")
print(f"{'-'*95}")

for cell in cell_types:
    raw_p = p_values_dict[cell]
    adj_q = corrected_p_dict[cell]
    
    # Logic to determine the change
    if adj_q < 0.05:
        status = "✅ SIGNIFICANT"
    elif raw_p < 0.05 and adj_q >= 0.05:
        status = "⚠️ LOST (Nominal Trend)"
    else:
        status = "❌ Not Significant"
        
    print(f"{cell:<20} | {raw_p:<12.4f} | {adj_q:<12.4f} | {status}")

print(f"{'='*95}")
