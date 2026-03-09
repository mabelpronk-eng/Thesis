#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# subset_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Subset for immune and myeloid cells
# *Adapated to create 3 subsets: myeloid, immune and lymphocyte 
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl) 
# Adapation: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# Usage:
"""
     python3 scripts/subset_scRNA.py \
     -i {input.adata_ingested_scaled} \
     -o1 {output.immune} \
     -o2 {output.myeloid} \
     -o3 {output.lymphocyte} \
     -cache_dir {params.cache_dir} 
"""

#================================================================================
# 0 Import libraries
#================================================================================
import argparse
import scanpy as sc

#================================================================================
# 1 Parse arguments
#================================================================================
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 subset_scRNA.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        'Subset for immune and myeloid cells')
    parser.add_argument('-i', help='path to input file',
                        dest='input',
                        type=str)
    parser.add_argument('-o1', help='path to immune subset output file',
                        dest='output1',
                        type=str)
    parser.add_argument('-o2', help='path to myeloid subset output file',
                        dest='output2',
                        type=str)
    parser.add_argument('-o3', help='path to lymphocyte subset output file',
                        dest='output3',
                        type=str)
    parser.add_argument('-cache_dir', help='Directory to save cache',
                        dest='cache_dir',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

# Set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

#===============================================================================
# 2 Subset
#===============================================================================
#-------------------------------------------------------------------------------
# 2.1 Read adata object
#-------------------------------------------------------------------------------
print("Reading " + args.input + "...")
adata = sc.read_h5ad(args.input)

#-------------------------------------------------------------------------------
# 2.2 Subset
#-------------------------------------------------------------------------------
# Subset for myeloid cells
myeloid_cells = adata[adata.obs['level1_celltype'].isin(['Myeloid', 'Microglia', 'TAM', 'Prolif_TAM','DC'])]
print(myeloid_cells)

# Subset for immune cells
immune_cells = adata[adata.obs['level0_celltype'].isin(['Immune'])]
print(immune_cells)

lymphocyte_cells = adata[adata.obs['level1_celltype'].isin(['B_cell', 'Plasma_B', 'T_NK', 'T_cell','T_reg','NK_cell'])]

#===============================================================================
# 3 Save subsets of adata
#===============================================================================
print("Writing adata immune...")
immune_cells.write(args.output1)

print("Writing adata myeloid...")
myeloid_cells.write(args.output2)

print("Writing adata lymphocytes...")
lymphocyte_cells.write(args.output3)
