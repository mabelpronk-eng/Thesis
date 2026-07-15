"""
Script: Export Statescope Signature Matrix

Description:
This script exports the cell-type signature matrix from a previously
initialized Statescope model.

The exported signature matrix contains the average expression and variance
for each gene across all cell types, together with an indicator specifying
whether a gene was selected as a marker by Statescope.

The script:
1. Loads a previously initialized Statescope model.
2. Extracts the cell-type expression and variance matrices.
3. Combines both matrices into a single signature table.
4. Annotates each gene as a marker or non-marker.
5. Saves the signature matrix as a tab-separated (.tsv) file for downstream
   analyses or inspection.

Inputs:
    --model    Initialized Statescope model (.pkl)
    --out      Output signature matrix (.tsv)

Output:
    Tab-separated signature matrix containing:
        - Gene name
        - Marker status
        - Mean expression per cell type
        - Expression variance per cell type

Example:
    python export_signature.py \
        --model statescope_model.pkl \
        --out signature_matrix.tsv

Adaptation: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""

import os
import sys
import pandas as pd
import argparse

# -------------------------------------------------
# Add the Statescope source directory to Python path
# -------------------------------------------------
# Assumes the following project structure:
# Project/
# ├── scripts/
# │   └── export_signature.py
# └── src/
#     └── Statescope/
THIS_DIR = os.path.dirname(__file__)
SRC_DIR = os.path.abspath(os.path.join(THIS_DIR, '..', 'src'))

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from Statescope.Statescope import Statescope

def main():

    # -------------------------------------------------
    # Parse command-line arguments
    # -------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    # -------------------------------------------------
    # Load the initialized Statescope model
    # -------------------------------------------------
    model = Statescope.load(args.model)

    # -------------------------------------------------
    # Construct the signature matrix
    # -------------------------------------------------
    # Combine mean expression and variance matrices
    signature = pd.concat([model.scExp, model.scVar], axis=1) 

    # Add gene names and indicate marker genes
    signature["Gene"] = signature.index 
    signature["IsMarker"] = signature["Gene"].isin(model.Markers) 

    # Reorder columns for readability
    cols = ["Gene", "IsMarker"] + [c for c in signature.columns if c not in ["Gene", "IsMarker"]] 
    signature = signature[cols] 

    # -------------------------------------------------
    # Save signature matrix
    # -------------------------------------------------
    signature.to_csv(args.out, sep="\t", index=False)
    print(f"Signature saved to {args.out}")
if __name__ == "__main__":
    main()
