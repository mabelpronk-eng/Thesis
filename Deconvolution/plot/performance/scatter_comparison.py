"""
Script: Statescope vs Ground Truth Scatterplot (Pseudobulk Evaluation)

Description:
This script evaluates deconvolution performance by comparing predicted
cell-type fractions (Statescope/CIBERSORTx output) against ground-truth
fractions derived from pseudobulk samples.

The script:
1. Loads predicted and ground-truth cell-type fraction matrices.
2. Aligns both datasets by shared samples and cell types.
3. Converts matrices to NumPy arrays for evaluation.
4. Computes performance metrics (Pearson correlation and RMSE).
5. Generates a per-cell-type scatterplot (predicted vs true fractions).
6. Saves a color mapping dictionary for consistent visualization across plots.
7. Outputs a publication-ready figure with overall performance metrics.

Outputs:
- Scatterplot (PDF) of predicted vs true cell fractions
- JSON file containing consistent cell-type color mapping
"""
# script based on scatter in script: /net/beegfs/users/P086608/StatescopePro_v2/tutorial/Output_pseudobulk/ELBO_decomposition_compare_models.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os 
print(os.getcwd())
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

# ==================================================
# 1. LOAD DATA
# ==================================================
# Predicted fractions (Statescope / CIBERSORTx)
#fractions_path = "/net/beegfs/users/P086608/Statescope/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_0001_rep_10/fractions3.csv"          # CIBERSORT output

# Ground-truth pseudobulk fractions
gt_path = "/net/beegfs/users/P086608/pseudobulk/level3_celltype/cell_type_fractions.csv"              # GT fractions

df_cibersort = pd.read_csv(fractions_path, sep = ',', index_col=0)
df_gt = pd.read_csv(gt_path, sep=",", index_col=0)

print("Loaded CIBERSORT:", df_cibersort.shape)
print("Loaded GT:", df_gt.shape)

# ==================================================
# 2. ALIGN SAMPLES + CELL TYPES
# ==================================================

# match cell types
common_celltypes = df_gt.columns.intersection(df_cibersort.columns)

df_cibersort = df_cibersort[common_celltypes]
df_gt = df_gt[common_celltypes]

# match samples
common_samples = df_gt.index.intersection(df_cibersort.index)

df_cibersort = df_cibersort.loc[common_samples]
df_gt = df_gt.loc[common_samples]

print("After alignment:")
print(df_cibersort.shape, df_gt.shape)

# ==================================================
# 3. CONVERT TO NUMPY / TENSORS
# ==================================================

X = df_cibersort.values   # predicted
Y = df_gt.values          # ground truth

cell_order = list(common_celltypes)
n_sample = X.shape[0]

# ==================================================
# 4. METRICS FUNCTIONS
# ==================================================

def safe_pcc(x, y):
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return np.corrcoef(x, y)[0, 1]

def safe_rmse(x, y):
    return np.sqrt(np.mean((np.asarray(x) - np.asarray(y)) ** 2))

# ==================================================
# 5. PLOTTING FUNCTION (same style as Statescope)
# ==================================================

def plot_scatter(x_np, y_np, cell_order, n_sample, title, out_file):

    x = x_np.ravel()
    y = y_np.ravel()

    overall_pcc = safe_pcc(x, y)
    overall_rmse = safe_rmse(x, y)

    ct_labels = np.tile(np.array(cell_order), n_sample)
    cmap = plt.get_cmap("tab20")
    ct_to_color = {ct: cmap(i % 20) for i, ct in enumerate(cell_order)}
    with open("celltype_colors.json", "w") as f:
        json.dump(ct_to_color, f, indent=4) #So i can use these colours later
    plt.figure(figsize=(7, 5))

    for j, ct in enumerate(cell_order):
        xm = x_np[:, j]
        ym = y_np[:, j]
        rmse_ct = safe_rmse(xm, ym)
        pcc_ct = safe_pcc(xm, ym)
        #m = (ct_labels == ct)

        plt.scatter(
            xm, ym,
            s=18, alpha=0.7,
            #label=f"{ct} (RMSE={rmse_ct:.3f}, PCC={pcc_ct:.3f})",
            label=ct,
            c=[ct_to_color[ct]]
        )

    lo = min(np.min(x), np.min(y))
    hi = max(np.max(x), np.max(y))
    plt.plot([lo, hi], [lo, hi], linewidth=1)

    plt.title(f"{title}\nPCC={overall_pcc:.3f}, RMSE={overall_rmse:.3f}")
    plt.xlabel("Predicted cell fractions")
    plt.ylabel("True cell fractions")
    plt.legend(title="Cell types", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    plt.xlim(0, 1)
    plt.ylim(0, 1)

    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()

# ==================================================
# 6. RUN PLOT
# ==================================================

plot_scatter(
    X,
    Y,
    cell_order,
    n_sample,
    title="Scatterplot of true vs predicted fractions",
    out_file="/net/beegfs/users/P086608/Statescope/StatescopePro_v2/tutorial/Output_pseudobulk/lam_0_001_rep_10/statescope_0.001_vs_gt.pdf"
)

print("Done → saved plot: cibersort_vs_gt.png")
