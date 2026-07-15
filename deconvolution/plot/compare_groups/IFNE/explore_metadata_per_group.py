"""
Clinical Covariate Analysis Pipeline

This script performs a systematic comparison of clinical characteristics across CNA-defined groups.
It integrates clinical metadata with group classifications and evaluates whether groups differ
in key covariates, which could indicate potential confounding effects.

Main steps:
1. Load and merge clinical and group data
2. Exclude samples with missing or incomplete clinical information
3. Summarize cohort composition per group
4. Assess distribution of categorical and continuous covariates
5. Visualize group differences
6. Perform statistical testing (focused on Group 3 vs Group 4)

The goal is to determine whether observed biological differences between groups
(e.g., TME composition) may be influenced by underlying clinical variables.

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- 1. LOAD DATA ---
# Load clinical metadata and CNA-based group assignments
df_meta = pd.read_csv('/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/clinical/GBM_Clinical_data.tsv', sep='\t')
df_groups = pd.read_csv('/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/final/classification_with_groups.csv')

# --- 2. PREP FOR JOIN ---
# Standardize patient/sample identifiers to ensure correct merging
# Remove leading/trailing whitespace in ID columns
df_meta['bcr_patient_barcode'] = df_meta['bcr_patient_barcode'].str.strip()
df_groups['SAMPLES'] = df_groups['SAMPLES'].str.strip()

# Samples to exclude due to missing or incomplete clinical data
# These are removed prior to merging to avoid bias in downstream analyses
exclude_samples = ['TCGA-06-0221', 'TCGA-14-0736', 'TCGA-14-1402', 'TCGA-19-0957']
df_groups = df_groups[~df_groups['SAMPLES'].isin(exclude_samples)]

# Merge clinical metadata with group assignments based on patient IDs
df_merged = pd.merge(df_groups, df_meta, left_on='SAMPLES', right_on='bcr_patient_barcode', how='inner')

# --- 2.5 SAMPLE COUNT SUMMARY ---
# Verify final sample size per group after filtering and merging
group_counts = df_merged['Group'].value_counts().sort_index()

print("\n" + "="*30)
print("FINAL SAMPLE COUNTS PER GROUP")
print("="*30)
print(group_counts)
print(f"Total Samples Matched: {len(df_merged)}")
print("="*30 + "\n")

# --- 3. COVARIATE AUDIT ---
# Evaluate distribution of categorical and continuous clinical variables across groups
# Define the variables of interest
categorical_covs = ['gender', 'prior_treatment', 'site_of_resection_or_biopsy', 'race','method_of_diagnosis']
continuous_covs = ['age_at_index']

# Generate normalized contingency tables (percentages per group)
print("--- Distribution of Covariates across Groups ---")

for cov in categorical_covs:
    print(f"\nSummary for {cov}:")
    # This creates a cross-tabulation table
    ctab = pd.crosstab(df_merged['Group'], df_merged[cov], normalize='index') * 100
    print(ctab.round(2))

# Summary statistics for age (continuous variable)
age_summary = df_merged.groupby('Group')['age_at_index'].agg(['mean', 'std', 'median'])
print("\nAge Summary per Group:")
print(age_summary)

# --- 4. VISUALIZATION ---
# Visualize age distribution across groups using boxplots
plt.figure(figsize=(10, 6))
sns.boxplot(data=df_merged, x='Group', y='age_at_index', palette='Set2')
plt.title('Age Distribution Across CNA Groups')
plt.ylabel('Age at Index')
save_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/clinical/boxplot_age_group.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight') 

# --- 5. PLOTTING FUNCTION FOR STACKED BAR CHARTS ---
# Generates stacked bar plots showing percentage distribution per group
def plot_stacked_bar(covariate, save_name):
    # Create percentage table
    ctab = pd.crosstab(df_merged['Group'], df_merged[covariate], normalize='index') * 100
    
    # Plot
    ax = ctab.plot(kind='bar', stacked=True, figsize=(10, 6))
    
    plt.title(f'{covariate.replace("_", " ").title()} Distribution Across Groups')
    plt.ylabel('Percentage (%)')
    plt.xlabel('Group')
    plt.legend(title=covariate, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f'/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/clinical/{save_name}', dpi=300)
    plt.close()

# --- 6. GENERATE PLOTS ---
plot_stacked_bar('race', 'race_distribution.png')
plot_stacked_bar('method_of_diagnosis', 'method_diagnosis_distribution.png')
plot_stacked_bar('gender', 'gender_distribution.png')

# Perform statistical tests to assess differences between Group 3 and Group 4
# Fisher’s Exact Test is used for 2x2 tables; Chi-square otherwise
# --- FILTER ONLY GROUP 3 AND 4 ---
df_34 = df_merged[df_merged['Group'].isin(['Group 3', 'Group 4'])]

from scipy.stats import chi2_contingency, fisher_exact

def test_group3_vs_4(covariate):
    print(f"\n--- {covariate.upper()} (Group 3 vs Group 4) ---")
    
    table = pd.crosstab(df_34['Group'], df_34[covariate])
    print("\nContingency Table:")
    print(table)

    # If 2x2 → use Fisher
    if table.shape == (2, 2):
        _, p = fisher_exact(table)
        print(f"Fisher’s Exact Test p-value = {p:.5f}")
    else:
        chi2, p, dof, expected = chi2_contingency(table)
        print(f"Chi-square p-value = {p:.5f}")
        
        if (expected < 5).sum() > 0:
            print("⚠️ Some expected counts < 5 → interpret cautiously")

# Run tests
test_group3_vs_4('gender')
test_group3_vs_4('race')
test_group3_vs_4('method_of_diagnosis')

# Test age differences between groups with assumption checking:
# - Normality (Shapiro-Wilk)
# - Homogeneity of variance (Levene)
# Select appropriate test (T-test or Mann-Whitney U) based on assumptions
from scipy.stats import shapiro, levene, mannwhitneyu, ttest_ind

def test_age_with_assumptions():
    age_g3 = df_34[df_34['Group'] == 'Group 3']['age_at_index'].dropna()
    age_g4 = df_34[df_34['Group'] == 'Group 4']['age_at_index'].dropna()

    # 1. Check Normality (Shapiro-Wilk)
    _, p_shap3 = shapiro(age_g3)
    _, p_shap4 = shapiro(age_g4)
    
    # 2. Check Variance (Levene)
    _, p_levene = levene(age_g3, age_g4)

    print(f"Shapiro P-values: Group3={p_shap3:.4f}, Group4={p_shap4:.4f}")
    print(f"Levene P-value: {p_levene:.4f}")

    # 3. Decision Logic
    if p_shap3 > 0.05 and p_shap4 > 0.05 and p_levene > 0.05:
        print("Data is Normal and Homogeneous -> Using T-Test")
        stat, p = ttest_ind(age_g3, age_g4)
    else:
        print("Assumptions failed -> Using Mann-Whitney U")
        stat, p = mannwhitneyu(age_g3, age_g4)
    
    print(f"Final Test P-value: {p:.5f}")
test_age_with_assumptions()
