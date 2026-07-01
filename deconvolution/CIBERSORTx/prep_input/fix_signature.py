# ============================================================
# Prepare single-cell signature of Statescope for CIBERSORT /
# deconvolution workflows
#
# Purpose:
# This script converts the generated single-cell signature
# made by Statescope back from log-normalized space to linear expression
# space so that it matches the preprocessing scale of the
# bulk RNA-seq input used for deconvolution.
#
# Original signature preprocessing:
#   sc.pp.normalize_total(target_sum=1e4)
#   sc.pp.log1p()
#
# Therefore, scExp values in the signature are:
#   log(CP10K + 1)
#
# This script:
#   1. Keeps only marker genes
#   2. Retains only scExp expression columns
#   3. Reverts log1p transformation using:
#
#         linear = exp(x) - 1
#
# Result:
# The final signature matrix is in linear CP10K space,
# matching the bulk RNA-seq preprocessing used for:
#
#   - CIBERSORT
#
# IMPORTANT:
# - Output remains in linear space
# - No additional scaling or log transformation should be applied
# - Bulk and signature matrices must remain in the same scale
#
# Author: Mabel Pronk

import pandas as pd
import numpy as np

# Load signature
sig = pd.read_csv("/net/beegfs/users/P086608/CIBERSORT/data/signature_17celltypes_GBM.txt", sep="\t")

# 1. Keep only marker genes
sig = sig[sig["IsMarker"] == True].copy()

# 2. Keep Gene + scExp columns
exp_cols = [c for c in sig.columns if c.startswith("scExp_")]
sig_linear = sig[["Gene"] + exp_cols].copy()

# 3. Convert log1p back to linear
sig_linear[exp_cols] = np.expm1(sig_linear[exp_cols])

# (optional) safety: replace tiny negatives from numerical error
#sig_linear[exp_cols] = sig_linear[exp_cols].clip(lower=0)

# 4. Save
sig_linear.to_csv("/net/beegfs/users/P086608/CIBERSORT/data/signature_CIBERSORTx_linear.tsv", sep="\t", index=False)