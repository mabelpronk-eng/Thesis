#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# normalize_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Normalize and log1p transform the data without scaling (due to ingestion)
# *Adapted to store filtered (QC) cells in an adata layer before transforming the data. 
#  Additionally a layer is added in which the log transformed and normalized data is stored (thus easier accessible later on)
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
# Adaptation: Mabel Pronk (m.pronk3@amsterdamumc.nl)
#

"""
        python3 scripts/normalize_scRNA.py \
        -i {input.adata_concatenated} \
        -o {output.adata_normalized} \
        -plot_dir {params.plot_dir} \
        -cache_dir {params.cache_dir}  
              
"""

#================================================================================================================
# 0 Import libraries
#================================================================================================================
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
import scrublet as scr

#================================================================================================================
# 1 Parse arguments
#================================================================================================================
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog='python3 normalize_scRNA.py',
                                     formatter_class=argparse.RawTextHelpFormatter, description='Normalization and log1p transformation of scRNA')
    parser.add_argument('-i', help='Input concatenated filtered adata', dest='input', type=str, required=True)
    parser.add_argument('-o', help='Path to output file', dest='output', type=str, required=True)
    parser.add_argument('-plot_dir', help='Directory to save plots', dest='plot_dir', type=str, required=True)
    parser.add_argument('-cache_dir', help='Directory to save cache', dest='cache_dir', type=str, required=True)
    args = parser.parse_args()
    return args

args = parse_args()
adata = sc.read_h5ad(args.input)

# Create a new layer called 'filtered' containing the current X matrix
adata.layers["filtered"] = adata.X.copy()

# Set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

#================================================================================================================
# 2 Preprocessing
#================================================================================================================
#-------------------------------------------------------------------------------
# 2.1 Normalize by library size and log transform
#-------------------------------------------------------------------------------
print("Normalizing by library size and log-transform")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store log-normalized values in a new layer instead of raw
adata.layers["log_norm"] = adata.X.copy()  # keep the log-normalized matrix

#-------------------------------------------------------------------------------
# 2.2 Identify highly variable genes
#-------------------------------------------------------------------------------
print("Identify and plot highly variable genes")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=4, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show = False)
plt.savefig(args.plot_dir + "total_counts_vs_pct_mt.png")

print(adata.var['highly_variable'].sum())

#-------------------------------------------------------------------------------
# 2.3 Select highly variable genes
#-------------------------------------------------------------------------------
adata = adata[:, adata.var.highly_variable]

#-------------------------------------------------------------------------------
# 2.4 Regress out effects of total_counts and pct_counts_mt (and scale to unit variance)
#-------------------------------------------------------------------------------
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#sc.pp.scale(adata, max_value=10) # leave out scaling due to subsequent ingestion

#================================================================================================================
# 3 Write to file
#================================================================================================================
adata.write(args.output)
