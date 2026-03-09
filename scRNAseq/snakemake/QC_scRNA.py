#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# QC_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# QC for scRNAseq counts of each dataset separately before Normalization
# *Adapted for glioma data*
#
# Author: Shiva Najjary (s.najjary@amsterdamumc.nl)
# Aadaptation: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl) & Mabel Pronk (m.pronk3@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/QC_scRNA.py \
        -i {input.matrix} \
        -o {output.adata_processed} \
        -min_genes {params.min_genes} \
        -max_genes {params.max_genes} \
        -min_cells {params.min_cells} \
        -max_mt_pct {params.max_mt_pct} \
        -cache_dir {params.cache_dir} \
        -plot_dir {params.plot_dir}
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
# 1 Parse arguments and prepare working directory
#================================================================================================================
#-------------------------------------------------------------------------------
# 1.1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog='python3 QC_scRNA.py',
                                     formatter_class=argparse.RawTextHelpFormatter, description='QC scRNAseq counts using scanpy')
    parser.add_argument('-i', help='Input matrix', dest='input', type=str, required=True)
    parser.add_argument('-o', help='Path to output file', dest='output', type=str, required=True)
    parser.add_argument('-min_genes', help='Minimal genes per cell', dest='min_genes', type=int, default=300)
    parser.add_argument('-max_genes', help='Maximal genes per cell', dest='max_genes', type=int, default=8000)
    parser.add_argument('-max_mt_pct', help='Maximal percentage of mitochondrial genes', dest='max_mt_pct', type=int, default=10)
    parser.add_argument('-min_cells', help='Minimal cells', dest='min_cells', type=int, default=50)
    parser.add_argument('-cache_dir', help='Directory to save cache', dest='cache_dir', type=str, required=True)
    parser.add_argument('-plot_dir', help='Directory to save plots', dest='plot_dir', type=str, required=True)
    args = parser.parse_args()
    return args

args = parse_args()

#-------------------------------------------------------------------------------
# 1.2 Configure scanpy and prepare working directory
#-------------------------------------------------------------------------------
# Set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

# Create plot directories
if not os.path.exists(args.plot_dir):
    os.makedirs(args.plot_dir)
if not os.path.exists(os.path.join(args.plot_dir, "before_filtering")):
    os.makedirs(os.path.join(args.plot_dir, "before_filtering"))
if not os.path.exists(os.path.join(args.plot_dir, "after_filtering")):
    os.makedirs(os.path.join(args.plot_dir, "after_filtering"))


#================================================================================================================
# 2 QC
#================================================================================================================
#----------------------------------------------------------------------------------------------------------------
# 2.1 Read matrix data to adata object
#----------------------------------------------------------------------------------------------------------------
# Read the 10x matrix data
print("Reading scRNAseq counts from directory")
adata = sc.read_h5ad(args.input)

print(adata)

# Make a copy for initial filtering
adata_copy = adata.copy()


#---------------------------------------------------------------------------------------------------------------------------------------------
# 2.2 Plot top 50 expressed genes and QC metrics
#---------------------------------------------------------------------------------------------------------------------------------------------
print("Plotting and saving top 50 expressed genes before normalization in " + args.plot_dir + "Top50_expressed_genes_prenorm.png")
sc.pl.highest_expr_genes(adata, n_top=50,show = False)
plt.savefig(args.plot_dir + "Top50_expressed_genes_prenorm.png")

#--------------------------------------------------------------------------------------------------------------------------------------------
# 2.3 Calculate/Plot QC metrics
#--------------------------------------------------------------------------------------------------------------------------------------------
# Annotate mitochrondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print("Plotting and saving QC metrics before filtering in "+ args.plot_dir + "before_filtering/")

# Plot number of genes, total counts, percentage mitochondrial genes
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show = False)
plt.savefig(args.plot_dir + "before_filtering/" + "QC_metrics.png")

# Plot percentage of mitochondrial counts vs total counts
print("Plotting percentage mitochondrial counts vs total counts")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show = False)
plt.savefig(args.plot_dir + "before_filtering/" + "total_counts_vs_pct_mt.png")

# Plot number of genes by counts vs total counts
print("Plotting and saving n_genes_by_counts vs total_counts")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show = False)
plt.savefig(args.plot_dir + "before_filtering/" + "total_counts_vs_ngenes.png")


#================================================================================================================
# 3 Filtering
#================================================================================================================
# Filter cells by min/max genes per cell
print(f"Filtering cells: min genes = {args.min_genes}, max genes = {args.max_genes}")
sc.pp.filter_cells(adata, min_genes=args.min_genes)
sc.pp.filter_cells(adata, max_genes=args.max_genes)
print(adata.n_obs, "cells remaining")

# Filter genes by min number of cells
print(f"Filtering genes: min cells = {args.min_cells}")
sc.pp.filter_genes(adata, min_cells=args.min_cells)
print(adata.n_obs, "cells remaining")


#================================================================================================================
# 4 Recalculate QC metrics and filter mitochondrial gene expression
#================================================================================================================
#----------------------------------------------------------------------------------------------------------------
# 4.1 Calculate mt_counts after filtering
#----------------------------------------------------------------------------------------------------------------
print("Recalculating QC metrics after filtering")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


#----------------------------------------------------------------------------------------------------------------
# 4.2 Filter cells by mitochondrial gene expression
#----------------------------------------------------------------------------------------------------------------
print("Retain cells with a maximum precentage of mitochondrial genes of",str(args.max_mt_pct) + "%")
adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :]
print(adata.n_obs, "cells remaining")

#----------------------------------------------------------------------------------------------------------------
# 4.3 Plot QC metrics after filtering
#----------------------------------------------------------------------------------------------------------------
print("Plotting and saving QC metrics after filtering in "+  args.plot_dir + "after_filtering/")
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show = False)
plt.savefig(args.plot_dir + "after_filtering/" + "QC_metrics.png")

# Plot percentage of mitochondrial counts vs total counts after filtering
print("Plotting percentage mitochondrial counts vs total counts after filtering")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show = False)
plt.savefig(args.plot_dir + "after_filtering/" + "total_counts_vs_pct_mt.png")

# Plot number of genes vs total counts after filtering
print("Plotting and saving n_genes_by_counts vs total_counts")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show = False)
plt.savefig(args.plot_dir + "after_filtering/" + "total_counts_vs_ngenes.png")


#================================================================================================================
# 5 Check data
#================================================================================================================
print("Filtered data summary:")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# Store raw counts
print("Storing filtered counts...")
adata.raw = adata

#================================================================================================================
# 6 Write filtered adata file
#================================================================================================================
print("Writing to file")
adata.write(args.output)
