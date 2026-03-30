"""
Script: Validation of Statescope-Derived Malignant Fractions Against Prior Purity Estimates

Description:
This script evaluates the concordance between malignant cell fractions estimated by 
Statescope (from bulk RNA-seq deconvolution) and independent prior purity estimates 
(e.g., derived from DNA copy number data). 

The script:
1. Loads deconvolution output and prior purity data.
2. Identifies mismatches in sample identifiers between datasets.
3. Restricts analysis to shared samples.
4. Computes Pearson correlation between the two purity estimates.
5. Generates a scatter plot with regression and identity (y=x) lines to assess agreement.

The resulting figure is saved for downstream reporting and validation purposes.

Author: Mabel Pronk

"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os

#-------------------------------------------------------------------------------
# 1. Load the Data
#-------------------------------------------------------------------------------

# Deconvolution fractions (Statescope output)
# Note: index_col=0 ensures sample IDs are used as row indices
deconv_path = '/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output/fractions3.csv'
df_deconv = pd.read_csv(deconv_path, index_col=0)

# Prior purity estimates (e.g., DNA-seq or pathology-based)
purity_path = '/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/input/malignant_fraction_TCGA.csv'
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
save_dir = "/net/beegfs/users/P086608/StatescopePro_v2/TCGA_bulk/Output"
save_path = os.path.join(save_dir, "purity_validation_correlation.png")

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
    
    #-------------------------------------------------------------------------------
    # 4. Visualization: Scatter Plot with Regression and Identity Line
    #-------------------------------------------------------------------------------
    
    plt.figure(figsize=(8, 7))
    
    # Plot data points with regression line
    sns.regplot(
        data=merged, 
        x='Malignant_Prior', 
        y='Malignant_Statescope',
        scatter_kws={'alpha':0.5, 's':40},
        line_kws={'color':'red', 'label': f'Regression (r={r:.2f})'}
    )
    
    # Add identity line (y = x) representing perfect agreement
    max_val = max(
        merged['Malignant_Prior'].max(), 
        merged['Malignant_Statescope'].max()
    )
    plt.plot([0, max_val], [0, max_val], color='black', linestyle='--', label='Identity (y=x)')
    
    # Customize plot appearance
    plt.title('Purity Validation: Trend vs. Accuracy')
    plt.xlabel('Prior Purity (DNA/Pathology)')
    plt.ylabel('Statescope Malignant Fraction')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    
    # Save figure
    plt.savefig(save_path, dpi=300)

else:
    print("\n[ERROR] No matching samples found. Check if the IDs (indices) format matches.")
