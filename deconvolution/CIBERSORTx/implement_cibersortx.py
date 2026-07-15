# ============================================================
# Script: CIBERSORT-style Deconvolution Using SVR + NNLS
# ============================================================
#
# Overview:
# This script implements a support vector regression (SVR)-based
# deconvolution approach inspired by CIBERSORT to estimate
# cell type proportions from bulk RNA-seq data using a reference
# single-cell-derived signature matrix.
#
# ------------------------------------------------------------
# Method Summary
# ------------------------------------------------------------
# The pipeline combines:
# 1. NuSVR (linear kernel) for robust regression-based estimation
#    of cell-type-specific expression contributions
# 2. Support vector selection (via minimum RMSE model selection)
# 3. Non-negative least squares (NNLS) refinement step
# 4. Post-hoc normalization to enforce compositional constraints
#
# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
# - Bulk expression matrix (genes × samples)
# - Signature matrix derived from scRNA-seq (genes × cell types)
#   (aligned on shared gene space)
#
# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
# - Cell type fraction matrix (samples × cell types)
# - Intermediate regression matrices (SVR coefficients, NNLS estimates)
# - Console validation of compositional constraints
#
#  Author: Yongsoo Kim (yo.kim@amsterdamumc.nl)
#  Adaptation: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# ============================================================

from sklearn.svm import SVR
from sklearn.svm import NuSVR
from sklearn.metrics import mean_squared_error as mse
import numpy as np
import pandas as pd
from scipy.optimize import nnls


def CIBERSORT(X, Y, Njob=1):
    Ngene, Nsample = Y.shape
    Ncell = X.shape[1]

    # estimate fraction
    SVRcoef = np.zeros((Ncell, Nsample))
    Selcoef = np.zeros((Ngene, Nsample))
    Nus = [0.25, 0.5, 0.75]
    for i in range(Nsample):
        sols = [NuSVR(kernel='linear', nu=nu).fit(X,Y[:,i]) for nu in Nus]
        RMSE = [mse(sol.predict(X), Y[:,i]) for sol in sols]
        Selcoef[sols[np.argmin(RMSE)].support_, i] = 1
        SVRcoef[:,i] = np.maximum(sols[np.argmin(RMSE)].coef_,0)

    # estimate per-cell expression
    NNLS_mat = np.zeros((Ngene, Ncell))
    for g in range(Ngene):
        NNLS_mat[g,:] = nnls(np.transpose(SVRcoef), Y[g,:])[0]

    return SVRcoef, NNLS_mat, Selcoef

# ==========================================
# 2. DATA LOADING & STRINGS ALIGNMENT
# ==========================================
print("Loading bulk data and reference signatures...")

# 1. Load data forcing the first column (column 0) to be the index
# This ignores whatever string name ('Gene', 'index', etc.) is at the top left corner
bulk_df = pd.read_csv("/net/beegfs/users/P086608/CIBERSORT/data/TCGA/input/bulk_cp10k.tsv", sep="\t", index_col=0)
sig_df = pd.read_csv("/net/beegfs/users/P086608/CIBERSORT/data/signature/signature_CIBERSORTx_linear.tsv", sep="\t", index_col=0)

print(f"Loaded Bulk Index Name: {bulk_df.index.name}, Shape: {bulk_df.shape}")
print(f"Loaded Sig Index Name: {sig_df.index.name}, Shape: {sig_df.shape}")

# Find intersecting genes present in both datasets
common_genes = bulk_df.index.intersection(sig_df.index)
print(f"Aligning matrices on {len(common_genes)} common marker genes.")

# Subset and sort to ensure strict matching row indexes
X_df = sig_df.loc[common_genes]
Y_df = bulk_df.loc[common_genes]

# Convert to raw numpy arrays for the scikit-learn function
X = X_df.values  # Shape: (Genes, CellTypes)
Y = Y_df.values  # Shape: (Genes, Samples)

# ==========================================
# 3. RUN THE PIPELINE
# ==========================================
SVRcoef, NNLS_mat, Selcoef = CIBERSORT(X, Y)
print("Deconvolution calculations finished!")

# ==========================================
# 4. FIX AFTER: CLIP NEGATIVES & NORMALIZE TO 1.0
# ==========================================
print("Applying post-processing corrections to regression coefficients...")

# Step 1: Convert to non-negative (negative values to zero)
SVRcoef_fixed = np.maximum(SVRcoef, 0)

# Step 2: Normalize to sum to one per sample (column-wise normalization)
# SVRcoef_fixed is shape (Ncell, Nsample), so we sum down the rows for each sample
sample_sums = SVRcoef_fixed.sum(axis=0)

# Avoid division by zero if a sample has all zeros
sample_sums[sample_sums == 0] = 1.0

# Divide each cell type's score by its sample total
SVRcoef_normalized = SVRcoef_fixed / sample_sums

# ==========================================
# 5. STORE AND FORMAT THE OUTPUT TABLES
# ==========================================
print("Formatting and saving output matrices...")


# Clean up the column names by removing "scExp_" prefix
clean_cell_types = [col.replace("scExp_", "") for col in X_df.columns]

# Create the clean fractions DataFrame
df_proportions = pd.DataFrame(SVRcoef_normalized, index=clean_cell_types, columns=Y_df.columns).T
df_proportions.index.name = "Sample_ID"

# --- SAFETY VALIDATION CHECK ---
# Print the row sums for the first 3 samples to prove they add up to 1.0 (100%)
row_sums = df_proportions.sum(axis=1)
print("\nValidation Check (Should all be 1.0):")
print(row_sums.head(3))

# Save the absolute true cell proportions file
output_path = "/net/beegfs/users/P086608/CIBERSORT/data/TCGA/output/cibersort_fractions_TCGAbulk_gbm.tsv"
df_proportions.to_csv(output_path, sep="\t")

print(f"\nSuccess! Clean cell fractions successfully saved to: {output_path}")
