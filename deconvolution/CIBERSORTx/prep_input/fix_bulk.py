#!/usr/bin/env python3
# ============================================================
# Create bulk RNA-seq input for CIBERSORT / deconvolution (works on pseudobulk and bulk from TCGA)
#
# Purpose:
# This script preprocesses bulk RNA-seq count data so that it
# matches the normalization scale used for the single-cell
# signature matrix.
#
# Signature preprocessing:
#   sc.pp.normalize_total(target_sum=1e4)
#   -> CP10K normalization
#
# Therefore, bulk RNA-seq data must also be transformed to:
#   Counts Per 10,000 (CP10K)
#
# IMPORTANT:
# - Input bulk data should be RAW COUNTS
# - Do NOT log-transform the bulk matrix
# - Output remains in linear CP10K space
#
# Output:
# A CP10K-normalized bulk expression matrix suitable for:
#   - CIBERSORT
#   - SVR-based deconvolution
#
# Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# ============================================================

import pandas as pd
import numpy as np

# ----------------------------
# 1. INPUT / OUTPUT PATHS
# ----------------------------
input_file  = "/net/beegfs/users/P086608/CIBERSORT/data/pseudobulk/pseudobulk_counts_transposed.csv"   # <-- change this
output_file = "/net/beegfs/users/P086608/CIBERSORT/data/pseudobulk/pseudobulk_cp10k.tsv"

# ----------------------------
# 2. LOAD DATA
# ----------------------------
# --- FIX 1: Change sep="\t" to sep="," (or removed it since comma is default for read_csv) if have bulk as csv file ---
bulk = pd.read_csv(input_file, sep=",")

# Ensure first column is gene names
bulk = bulk.set_index(bulk.columns[0])

# Convert everything to numeric (important!)
bulk = bulk.apply(pd.to_numeric, errors="coerce").fillna(0)

print(f"Loaded bulk matrix: {bulk.shape[0]} genes × {bulk.shape[1]} samples")

# ----------------------------
# 3. CP10K NORMALIZATION
# ----------------------------
# library-size normalize per sample (column-wise)
bulk_cp10k = bulk.div(bulk.sum(axis=0), axis=1) * 1e4

# safety check
print("Max per-sample sum after normalization:",
      bulk_cp10k.sum(axis=0).head())

# ----------------------------
# 4. SAVE OUTPUT
# ----------------------------
# --- FIX 2: Change sep="\t" to sep="," to save it back out as a true clean CSV file if wanted ---
bulk_cp10k.reset_index().to_csv(output_file, sep="\t", index=False)

print(f"Saved CP10K-normalized bulk matrix to: {output_file}")


