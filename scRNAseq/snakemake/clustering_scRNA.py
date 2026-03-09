#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# clustering_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Clustering and dimensionality reduction of scRNAseq counts with bbknn batch correction (possible to leave out)
# *Adapted to include possibility to perform Harmonay and BBKNN
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
# Adaptation: Mabel Pronk (m.pronk3@amsterdam.umc.nl)
#
# Usage:
"""
     python3 scripts/Clustering_scRNAseq_data.py \
Arguments updated:
    -i            Path to input AnnData (.h5ad) file
    -o            Path to output AnnData (.h5ad) file
    -n_pcs        Number of principal components for PCA (default: 25)
    -n_neighbors  Number of neighbors for KNN graph (default: 20)
    -resolution   Resolution parameter for Leiden clustering (default: 1)
    -cache_dir    Directory to store cached computations
    -plot_dir     Directory to save plots
    --bbknn       Optional flag to use BBKNN batch correction
    --harmony     Optional flag to use Harmony batch correction
"""

#================================================================================
# 0 Import libraries
#================================================================================
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import bbknn

#================================================================================
# 1 Parse arguments
#================================================================================
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 clustering_scRNAseq.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Perform clustering and visualization of scRNAseq data  ')
    parser.add_argument('-i', help='path to input file',
                        dest='input',
                        type=str)
    parser.add_argument('-o', help='path to output file',
                        dest='output',
                        type=str)
    parser.add_argument('-n_pcs', help='Number of Principial components',
                        dest='n_pcs',
                        default = 25,
                        type=int)
    parser.add_argument('-n_neighbors', help='Number of neighbors for KNN graph',
                        dest='n_neighbors',
                        default = 20,
                        type=int)
    parser.add_argument('-resolution', help='Clustering resolution',
                        dest='resolution',
                        default = 1,
                        type=float)
    parser.add_argument('-cache_dir', help='Directory to save cache',
                        dest='cache_dir',
                        type=str)
    parser.add_argument('-plot_dir', help='Directory to save plots',
                        dest='plot_dir',
                        type=str)
    parser.add_argument('--bbknn', help='Use BBKNN batch correction',
                    action='store_true')
    parser.add_argument('--harmony', help='Use Harmony batch correction',
                    action='store_true')

    args = parser.parse_args()
    return args

args = parse_args()

# Set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

if not os.path.exists(args.plot_dir): os.makedirs(args.plot_dir)

#===============================================================================
# 2 Clustering
#===============================================================================
#-------------------------------------------------------------------------------
# 2.1 Read adata object
#-------------------------------------------------------------------------------
print("Reading " + args.input + "...")
adata = sc.read_h5ad(args.input)

#-------------------------------------------------------------------------------
# 2.2 Perform PCA
#-------------------------------------------------------------------------------
sc.tl.pca(adata, svd_solver='arpack')

# Plot variance ratios of PCs
sc.pl.pca_variance_ratio(adata, log=False, n_pcs = 30, show= False)
plt.savefig(args.plot_dir + "PCA_variance_explained.png")


#===============================================================================
# 2.3 Compute neighbors (with or without batch correction)
#===============================================================================
print("Computing neighborhood graph...")

# --- BBKNN ---
if args.bbknn:
    print("→ Running BBKNN batch correction")
    sc.external.pp.bbknn(adata, batch_key="Cohort")

# --- Harmony ---
elif args.harmony:
    print("→ Running Harmony batch correction")

    # Run Harmony
    sc.external.pp.harmony_integrate(adata, key="Cohort")

    # Use Harmony-corrected PCs
    print("→ Using Harmony-corrected PCs for neighbors")
    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        use_rep="X_pca_harmony")

    

# --- No batch correction ---
else:
    print("→ No batch correction (standard neighbors)")
    print(f"   Using n_neighbors = {args.n_neighbors}, n_pcs = {args.n_pcs}")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
      



#-------------------------------------------------------------------------------
# 2.4 Compute UMAP space
#-------------------------------------------------------------------------------
print("Computing UMAP...")
sc.tl.umap(adata)

#-------------------------------------------------------------------------------
# 2.5 Perform Leiden clustering
#-------------------------------------------------------------------------------
print("Leiden clustering...")
print(f'   Using resolution = {args.resolution}')
sc.tl.leiden(adata, resolution = args.resolution)

#-------------------------------------------------------------------------------
# 2.6 Compute PAGA initialization
#-------------------------------------------------------------------------------
print("Computing PAGA...")
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)

#-------------------------------------------------------------------------------
# 2.7 Re-compute UMAP space with paga intialization
#-------------------------------------------------------------------------------
print("Recomputing UMAP...")
sc.tl.umap(adata, init_pos = 'paga')

#-------------------------------------------------------------------------------
# 2.8 Perform Leiden clustering on PAGA UMAP space
#-------------------------------------------------------------------------------
print("Leiden clustering...")
sc.tl.leiden(adata, resolution = args.resolution)


#===============================================================================
# 3 Save adata to file
#===============================================================================
print("Writing file...")
adata.write(args.output)
