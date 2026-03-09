#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ingest_annotations_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Transfer labels from the annotated data to the unannotated data
# *Adapted to work with different adata layer format 
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
# Adaptation: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# Usage:
"""
        python3 scripts/ingest_annotations_scRNA.py \
        -i {input.adata_normalized} \
        -o {output.adata_un_ingested_scaled} \
        -o_comb {output.adata_ingested_scaled} \
        -cache_dir {params.cache_dir}  
              
"""

#================================================================================================================
# 0 Import libraries
#================================================================================================================
import scanpy as sc
import numpy as np
import argparse

#================================================================================================================
# 1 Parse arguments
#================================================================================================================
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog='python3 ingest_annotations_scRNA.py',
                                     formatter_class=argparse.RawTextHelpFormatter, description='Transfer cell labels from annotated cells to unannotated cells')
    parser.add_argument('-i', help='Input normalized adata', dest='input', type=str, required=True)
    parser.add_argument('-o', help='Output file with transferred labels with only unannotated cells', dest='output', type=str, required=True)
    parser.add_argument('-o_comb', help='Output annotated ingested file with combined annotated and unannotated cells', dest='output_combined', type=str, required=True)
    parser.add_argument('-cache_dir', help='Directory to save cache', dest='cache_dir', type=str, required=True)
    args = parser.parse_args()
    return args

args = parse_args()

# Set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

#================================================================================================================
# 2 Ingest labels
#================================================================================================================
#-------------------------------------------------------------------------------
# 2.1 Split the data to annotated and unannotated cells
#-------------------------------------------------------------------------------
# Load the concatenated normalized adata
print("Loading adata...")
adata = sc.read_h5ad(args.input)

print(adata)
print(adata.obs)

# Split adata into reference and target (query) dataset based on the presence/absence of cell type label
print("Splitting adata...")
adata_ref = adata[~adata.obs['level1_celltype'].isna(), :]
adata_target = adata[adata.obs['level1_celltype'].isna(), :]

print(adata_ref)
print(adata_target)

#-------------------------------------------------------------------------------
# 2.2 Calculate PCA, neighbors and UMAP
#-------------------------------------------------------------------------------
print("Calculating PCA...")
sc.pp.pca(adata_ref, svd_solver='arpack')

print("Calculating neighbors...")
sc.pp.neighbors(adata_ref)

print("Computing UMAP...")
sc.tl.umap(adata_ref)

#-------------------------------------------------------------------------------
# 2.3 Transfer labels using the reference PCA space and kNN neighbors
#-------------------------------------------------------------------------------
print('Using ingest to map the target data...')
sc.tl.ingest(adata_target, adata_ref, obs=["level1_celltype", "level0_celltype"])
print(adata_target)

#-------------------------------------------------------------------------------
# 2.4 Save annotated target adata
#-------------------------------------------------------------------------------
print('Saving adata_target...')
adata_target.write(args.output)
print(adata_target)

#================================================================================================================
# 3 Concatenate and save adata to file
#================================================================================================================
#-------------------------------------------------------------------------------
# 3.1 Concatenate and scale
#-------------------------------------------------------------------------------
adata_concat = adata_ref.concatenate(adata_target, batch_categories=["ref", "target"])
print(adata_concat)
print(adata_concat.obs)

sc.pp.scale(adata_concat, max_value=10)

#-------------------------------------------------------------------------------
# 3.2 Save concatenated reference and target adata
#-------------------------------------------------------------------------------
adata_concat.write(args.output_combined)
