"""
Script: Initialize Statescope Model

Description:
This script initializes a Statescope deconvolution model using a bulk RNA-seq
expression matrix and a single-cell reference atlas.

The script:
1. Loads bulk RNA-seq expression data.
2. Loads the annotated single-cell reference atlas (.h5ad).
3. Initializes the Statescope model using the specified cell-type annotation.
4. Automatically derives marker genes and cell types from the reference atlas.
5. Prints the detected marker genes and cell types for quality control.
6. Saves the initialized model as a serialized (.pkl) object for downstream
   deconvolution.

Inputs:
    --bulk           Bulk RNA-seq expression matrix (.csv)
    --atlas          Annotated single-cell atlas (.h5ad)
    --celltype_key   Column in adata.obs containing cell-type annotations
    --out            Output filename (.pkl)

Output:
    Initialized Statescope model ready for deconvolution.

Example:
    python initialize.py \
        --bulk bulk_expression.csv \
        --atlas atlas.h5ad \
        --celltype_key lv1_celltype \
        --out statescope_model.pkl

Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""
import os
import sys
import argparse
import pandas as pd
import anndata

# -------------------------------------------------
# Add the Statescope source directory to Python path
# -------------------------------------------------
# Assumes the following project structure:
# Project/
# ├── scripts/
# │   └── initialize.py
# └── src/
#     └── Statescope/
# This assumes your structure is: Project/scripts/initialize.py 
# and the source is at: Project/src/Statescope

THIS_DIR = os.path.dirname(__file__)
SRC_DIR = os.path.abspath(os.path.join(THIS_DIR, '..', 'src'))

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from Statescope.Statescope import Initialize_Statescope

def main():
    # -------------------------------------------------
    # Parse command-line arguments
    # -------------------------------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument('--bulk', required=True)
    parser.add_argument('--atlas', required=True)
    parser.add_argument('--celltype_key', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    # -------------------------------------------------
    # Load bulk RNA-seq data and single-cell atlas
    # -------------------------------------------------
    Bulk = pd.read_csv(args.bulk, index_col=0)
    scRNAseq_dataset = anndata.read_h5ad(args.atlas)

    # -------------------------------------------------
    # Initialize the Statescope model
    # -------------------------------------------------
    model = Initialize_Statescope(Bulk, Signature=scRNAseq_dataset, celltype_key=args.celltype_key)

    # -------------------------------------------------
    # Print marker genes and detected cell types
    # -------------------------------------------------
    print('First 10 markers:',model.Markers[0:10]) # Extra genes can also be manually added: model.Markers = model.Markers + ['MARKERNAME1','MARKERNAME2']
    print('Celltypes:',model.Celltypes)

    # -------------------------------------------------
    # Save initialized model
    # -------------------------------------------------
    save_path = args.out.replace('.pkl', '')
    model.save(save_path, to_cpu=True)

if __name__ == "__main__":
    main()
