"""
===========================================================================
Cell Type Fraction and Malignant Expectation Calculator
===========================================================================

This script calculates per-sample cell type fractions from a single-cell 
AnnData (.h5ad) object and generates a CSV with expected fractions for 
deconvolution purposes.

Functionality:
1. Computes the fraction of each cell type per patient/sample.
2. Saves a CSV containing all cell type fractions.
3. Creates an "Expectation" CSV where only the 'Malignant' cell type is 
   kept; all other cell types are set to NaN.
4. Adjusts values of 0 and 1 for the 'Malignant' fraction to 0.01 and 
   0.99, respectively, to avoid issues in downstream deconvolution.

Inputs (command-line arguments):
1. Input .h5ad file containing single-cell expression data.
2. Column name in `adata.obs` for patient/sample IDs.
3. Column name in `adata.obs` for cell type labels.
4. Output CSV path for full cell type fractions.
5. Output CSV path for malignant fraction expectation.

Outputs:
- CSV file with all cell type fractions per sample.
- CSV file with malignant fraction expectation for deconvolution.

Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
Adaptation: Mabel Pronk (m.pronk3@amsterdamumc.nl)
"""

import sys
import pandas as pd
import numpy as np
import scanpy as sc



def calculate_fractions_expectation_csv(adata, patient_id, cell_type_column, cell_fractions_file, malignant_fraction_file):
    # Count the number of cells per cell type in each sample
    print("Calculating cell type fractions...")
    cell_type_counts = adata.obs.groupby([patient_id, cell_type_column]).size().unstack(fill_value=0)

    total_cells_per_sample = adata.obs.groupby(patient_id).size()   # Calculate the total number of cells per sample
    cell_type_fractions = cell_type_counts.div(total_cells_per_sample, axis=0)  # Compute the fraction of each cell type in each sample

    # Reorder the columns based on order of cell types in signature
    signature_order =  ['Malignant', 'TAM', 'Oligodendrocyte', 'Endothelial', 'Pericyte', 'CD4_Tcell', 'NK_cell', 'CD8_Tcell', 'Neutrophil', 'cDC', 'T_reg', 'Monocyte', 'B_cells', 'pDC', 'Plasma_B', 'Fibroblast']
    cell_type_fractions = cell_type_fractions[signature_order]
    print(cell_type_fractions.head())   # Print the first few rows of cell type fractions

    # Save the dataframe with all the cell type fractions
    print('Saving cell type fractions to CSV...')
    cell_type_fractions.to_csv(cell_fractions_file)

    # Keep values only for Malignant celltype, replace all other values with nan
    print("Calculating Expectation...")
    Expectation=cell_type_fractions
    Expectation.loc[:, Expectation.columns != 'Malignant'] = np.nan
    print(Expectation)

    Expectation = Expectation.map(lambda x: 0.01 if x == 0 else (0.99 if x == 1 else x), na_action='ignore')    # Replace 0 with 0.1 and 1 with 0.99 (necessary for deconvolution)

    # Save the dataframe with only malignant fraction and others nan to .csv
    print('Saving Expectations to CSV...')
    Expectation.to_csv(malignant_fraction_file)

    return cell_type_fractions, Expectation

if __name__ == "__main__":
    adata_file = sys.argv[1]                  # e.g., /net/beegfs/users/P086608/adata.h5ad
    patient_id_col = sys.argv[2]              # e.g., "patient_id"
    cell_type_col = sys.argv[3]               # e.g., "cell_type"
    cell_fractions_csv = sys.argv[4]          # e.g., "cell_type_fractions.csv"
    malignant_csv = sys.argv[5]               # e.g., "malignant_fraction.csv"

    adata = sc.read_h5ad(adata_file)

    calculate_fractions_expectation_csv(
        adata, patient_id_col, cell_type_col, cell_fractions_csv, malignant_csv)
