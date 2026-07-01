"""
====================================================================
Pseudobulk Absolute Error Comparison Between Statescope Models
====================================================================

Description:
------------
This script compares cell type fraction prediction performance between
two Statescope model versions using pseudobulk benchmark data.

For each cell type, the script:
    1. Loads predicted cell fractions from each model
    2. Compares predictions against ground truth fractions
    3. Computes absolute prediction errors
    4. Performs statistical testing between models
    5. Applies Benjamini–Hochberg FDR correction
    6. Visualizes error distributions with significance annotations

Statistical workflow:
---------------------
- Shapiro–Wilk test:
    Assess normality of each model distribution

- Levene’s test:
    Assess equality of variances

- Statistical test selection:
    * Student’s t-test        -> normal + equal variance
    * Welch’s t-test          -> normal + unequal variance
    * Mann–Whitney U test     -> non-normal distributions

- Multiple testing correction:
    Benjamini–Hochberg FDR correction across cell types

Outputs:
--------
1. Console summary of statistical test results
2. Annotated boxplot comparing absolute errors
3. Saved figure:
   Absolute_Error_Significance_Comparison.png

Author:
-------
Mabel Pronk

====================================================================
"""

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from scipy.stats import shapiro, levene, mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests
from statannotations.Annotator import Annotator

# --------------------------------------------------
# 1) Setup environment and import Statescope
# --------------------------------------------------
SRC_DIR = "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/src"
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from Statescope.Statescope import Statescope

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

# --------------------------------------------------
# 2) Define model and data paths
# --------------------------------------------------
model_paths = {
    "Original": "/net/beegfs/users/P086608/Statescope/StatescopePro_original/tutorial/pseudobulk_level3/Output_alpha_rem/statescope.pkl",
    "Optimized (Lambda 0.0001)": "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_0001_rep_10/statescope.pkl"
}

# Ground-truth cell type fractions
gt_path = "/net/beegfs/users/P086608/pseudobulk/level3_celltype/cell_type_fractions.csv"

# --------------------------------------------------
# 3) Load data and compute absolute errors
# --------------------------------------------------
# Load ground-truth fractions
gt_df = pd.read_csv(gt_path, index_col=0)

# Container for all per-sample errors
error_data = []

# Will store consistent cell type ordering
cell_order = []

# Loop through all model versions
for label, path in model_paths.items():

    # Load trained Statescope model
    model = Statescope.load(path, device="cpu")

    # Store cell type order from first model
    if not cell_order:
        cell_order = list(model.Celltypes)
    
    # Align GT columns to model cell order
    current_gt = gt_df.reindex(columns=cell_order).to_numpy()

    # Predicted cell fractions
    f_pred = model.BLADE.ExpF(model.BLADE.Beta).detach().cpu().numpy()

    # Compute absolute prediction error
    abs_error = np.abs(f_pred - current_gt)
    
    # Store error values per sample and cell type
    for i, ct in enumerate(cell_order):
        sample_errors = abs_error[:, i]
        for err in sample_errors:
            error_data.append({
                "Model": label,
                "Cell Type": ct,
                "Absolute Error": err
            })
# Convert to DataFrame for plotting/statistics
df_error = pd.DataFrame(error_data)

# --------------------------------------------------
# 4) Statistical testing per cell type
# --------------------------------------------------
raw_p_values = []

# Model labels
model_names = list(model_paths.keys()) # ["Original", "Optimized (Lambda 0.0001)"]

print(f"\n{'Cell Type':<20} | {'Test Used':<15} | {'Raw P-value'}")
print("-" * 50)

# Compare models independently for each cell type
for ct in cell_order:

    # Extract errors for both models
    data1 = df_error[(df_error['Cell Type'] == ct) & (df_error['Model'] == model_names[0])]['Absolute Error'].values
    data2 = df_error[(df_error['Cell Type'] == ct) & (df_error['Model'] == model_names[1])]['Absolute Error'].values
    
    # ----------------------------------------------
    # 4.1 Test normality (Shapiro–Wilk)
    # ----------------------------------------------
    _, p_norm1 = shapiro(data1)
    _, p_norm2 = shapiro(data2)
    is_normal = (p_norm1 > 0.05 and p_norm2 > 0.05)
    
    # ----------------------------------------------
    # 4.2 Test equal variances (Levene’s test)
    # ----------------------------------------------
    _, p_lev = levene(data1, data2)
    is_equal_var = (p_lev > 0.05)
    
    # ----------------------------------------------
    # 4.3 Select appropriate statistical test
    # ----------------------------------------------
    if is_normal and is_equal_var:
        test_name = "T-test"
        _, p_val = ttest_ind(data1, data2, equal_var=True)
    elif is_normal and not is_equal_var:
        test_name = "Welch's T"
        _, p_val = ttest_ind(data1, data2, equal_var=False)
    else:
        test_name = "Mann-Whitney"
        _, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
    
    # Store raw p-value
    raw_p_values.append(p_val)

    # Print summary
    print(f"{ct:<20} | {test_name:<15} | {p_val:.4e}")

# --------------------------------------------------
# 5) Multiple testing correction
# --------------------------------------------------
# Benjamini–Hochberg FDR correction
_, adj_p_values, _, _ = multipletests(raw_p_values, method='fdr_bh')

# --------------------------------------------------
# 6) Visualization
# --------------------------------------------------
plt.figure(figsize=(12, 5))
sns.set_style("ticks")

ax = sns.boxplot(
    data=df_error,
    x="Cell Type",
    y="Absolute Error",
    hue="Model",
    hue_order=model_names,
    palette=["#ca4370", "#3ca3e7"],
    linewidth=1.2,
    showfliers=False
)

# --------------------------------------------------
# 6.1 Add statistical significance annotations
# --------------------------------------------------
# Construct pairs for each cell type
pairs = [((ct, model_names[0]), (ct, model_names[1])) for ct in cell_order]

annotator = Annotator(ax, pairs, data=df_error, x="Cell Type", y="Absolute Error", hue="Model", hue_order=model_names)
annotator.configure(text_format="star", loc="inside")

# Use FDR-adjusted p-values
annotator.set_pvalues(adj_p_values) # Use the adjusted p-values
annotator.annotate()

# --------------------------------------------------
# 6.2 Plot styling
# --------------------------------------------------
#for spine in ax.spines.values():
#    spine.set_visible(True)
#    spine.set_color('black')
#    spine.set_linewidth(1.5)

#ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.3)

plt.title("Statistical Comparison of Absolute Error: Original vs. Optimized", pad=25)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.ylabel("Absolute Error")
plt.xlabel("Cell Types")
plt.xticks(rotation=45, ha='right')
plt.legend(title="Model Version", loc='upper right')

# Increase ylim slightly to prevent stars from being cut off
plt.ylim(0, df_error["Absolute Error"].max() * 1.2)

plt.tight_layout()

# --------------------------------------------------
# 7) Save figure
# --------------------------------------------------
output_fig = "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/tutorial/Output_pseudobulk/Absolute_Error_Significance_Comparison.png"
plt.savefig(output_fig, dpi=300)
print(f"\nAnnotated plot saved to: {output_fig}")
plt.show()
