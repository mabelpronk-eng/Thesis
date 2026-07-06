"""
Script: Validation of Deconvolution-Derived Malignant Fractions Against Prior Purity Estimates

Description:
This script evaluates the concordance between malignant cell fractions estimated by 
a deconvolution method (Statescope or CIBERSORTx) and independent prior purity estimates 
(e.g., derived from DNA copy number data such as ACE).

The script performs the following steps:
1. Loads deconvolution output and prior purity data.
2. Identifies mismatches in sample identifiers between datasets.
3. Restricts the analysis to shared samples.
4. Computes validation metrics:
   - Pearson correlation coefficient (PCC)
   - P-value
   - Root mean square deviation (RMSD)
5. Generates a scatter plot comparing predicted vs. prior tumor fractions,
   including an identity line (y = x) to assess agreement.
6. Compares correlation results obtained from SciPy and a custom NumPy implementation
   to ensure numerical consistency and robustness.

The resulting figure and metrics are used for downstream reporting and validation 
of deconvolution performance.

Author: Mabel Pronk
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os

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

#-------------------------------------------------------------------------------
# 1. Load the Data
#-------------------------------------------------------------------------------

# Deconvolution fractions (Statescope output)
# Note: index_col=0 ensures sample IDs are used as row indices
deconv_path = '/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output/statescope/fractions3.csv'
#deconv_path = '/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output/cibersort_fractions_TCGAbulk_gbm.tsv'
df_deconv = pd.read_csv(deconv_path, index_col=0) #if use tsv, add sep = '\t'

# Prior purity estimates (e.g., DNA-seq or pathology-based)
purity_path = '/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/input/malignant_fraction_TCGA.csv'
df_purity = pd.read_csv(purity_path, index_col=0)

#-------------------------------------------------------------------------------
# 2. Identify Non-Matching Samples
#-------------------------------------------------------------------------------

# Extract sample identifiers
deconv_idx = set(df_deconv.index)
purity_idx = set(df_purity.index)

# Determine overlap and mismatches
only_in_deconv = deconv_idx - purity_idx
only_in_purity = purity_idx - deconv_idx
shared_samples = deconv_idx.intersection(purity_idx)

print("--- Data Alignment Check ---")
print(f"Samples in Deconvolution: {len(deconv_idx)}")
print(f"Samples in Purity File:    {len(purity_idx)}")
print(f"Shared Samples (Overlap): {len(shared_samples)}")

# Report mismatched samples for troubleshooting
if only_in_deconv:
    print(f"\n[!] {len(only_in_deconv)} samples found in Deconvolution results but NOT in Purity file:")
    print(sorted(list(only_in_deconv)))

if only_in_purity:
    print(f"\n[!] {len(only_in_purity)} samples found in Purity file but NOT in Deconvolution results:")
    print(sorted(list(only_in_purity)))

#-------------------------------------------------------------------------------
# 3. Merge Data and Perform Correlation Analysis (Shared Samples Only)
#-------------------------------------------------------------------------------

# Define output path for the validation plot
save_dir = "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/TCGA_bulk/Output"
save_path = os.path.join(save_dir, "purity_validation_correlation.pdf")

if len(shared_samples) > 0:
    # Merge datasets on shared sample IDs
    # Assumes the relevant column in both datasets is labeled 'Malignant'
    merged = df_deconv[['Malignant']].merge(
        df_purity[['Malignant']], 
        left_index=True, 
        right_index=True, 
        suffixes=('_Statescope', '_Prior')
    ).dropna()

    # Compute Pearson correlation between estimates
    r, p = pearsonr(merged['Malignant_Statescope'], merged['Malignant_Prior'])
    rmsd = np.sqrt(((merged['Malignant_Statescope'] - merged['Malignant_Prior'])
                    
    #-------------------------------------------------------------------------------
    # 4. Visualization: Scatter Plot with Regression and Identity Line
    #-------------------------------------------------------------------------------
    
    plt.figure(figsize=(8, 7))
    
    # Plot data points with regression line
    sns.regplot(
        data=merged, 
        x='Malignant_Prior', 
        y='Malignant_Statescope',
        fit_reg=False,
        scatter_kws={'alpha':0.5, 's':40}
        #line_kws={'color':'red', 'label': f'Regression (r={r:.2f})'}
    )
    
    # Add identity line (y = x) representing perfect agreement
    max_val = max(
        merged['Malignant_Prior'].max(), 
        merged['Malignant_Statescope'].max()
    )
    plt.plot([0, max_val], [0, max_val], color='black', linestyle='--', label='Identity (y=x)')
    
    # Customize plot appearance
    plt.title('Scatterplot of true vs predicted tumor fraction\n'
    f'PCC = {r:.2f}, RMSD = {rmsd:.2f}')
    plt.xlabel('Prior tumor fraction (ACE)')
    plt.ylabel('Statescope tumor fraction estimate')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    
    # Save figure
    plt.savefig(save_path, dpi=300)

else:
    print("\n[ERROR] No matching samples found. Check if the IDs (indices) format matches.")


#-------------------------------------------------------------------------------
# 5. Quantitative Validation Metrics
#-------------------------------------------------------------------------------

import numpy as np 

# 1. Calculate Pearson Correlation (PCC)
# pearsonr returns (correlation coefficient, p-value)
r, p_val = pearsonr(merged['Malignant_Statescope'], merged['Malignant_Prior'])

# 2. Calculate RMSD (Root Mean Square Deviation)
# Measures the average deviation from perfect agreement (identity line y = x)
rmsd = np.sqrt(((merged['Malignant_Statescope'] - merged['Malignant_Prior']) ** 2).mean())

print(f"--- Validation Metrics (Malignant Cells) ---")
print(f"Pearson Correlation (r): {r:.4f}")
print(f"P-value:                {p_val:.4e}")
print(f"RMSD:                   {rmsd:.4f}")


#-------------------------------------------------------------------------------
# 6. Consistency Check: Correlation Implementation
#-------------------------------------------------------------------------------

def safe_pcc(x, y):
    """
    Computes Pearson correlation using NumPy with safeguards against zero variance.
    Returns NaN if either input has zero standard deviation.
    """
    x = np.asarray(x).reshape(-1)
    y = np.asarray(y).reshape(-1)
    
    # Prevent division by zero when variance is zero
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
        return np.nan
    
    return np.corrcoef(x, y)[0, 1]


# Calculate using SciPy (reference implementation)
r_scipy, p_val = pearsonr(merged['Malignant_Statescope'], merged['Malignant_Prior'])

# Calculate using custom NumPy implementation
r_np = safe_pcc(merged['Malignant_Statescope'], merged['Malignant_Prior'])

print(f"--- Correlation Comparison ---")
print(f"Scipy Pearson r: {r_scipy:.6f}")
print(f"NumPy safe_pcc:  {r_np:.6f}")
print(f"Difference:      {abs(r_scipy - r_np):.2e}")
