#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# label_scRNA.py Heslo@100
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Refine cell labels based on marker gene expression in immune and myeloid subset and in combined adata file
# 
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl) 
#
# Usage:
"""
     python3 scripts/label_scRNA.py \
     -i {input.adata_clustered} \
     -l {input.adata_labels} \
     -o {output.adata_labeled} \
     -cache_dir {params.cache_dir} 
"""

#============================================================================
# 0 Import libraries
#============================================================================
import argparse
import scanpy as sc
import pandas as pd

#============================================================================
# 1 Parse arguments
#============================================================================
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 label_scRNA.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        'Label cells in immune and myeloid subset')
    parser.add_argument('-i', help='path to input subset file',
                        dest='input',
                        type=str)
    parser.add_argument('-l', help='path to labels matched to leiden clusters',
                        dest='labels',
                        type=str)
    parser.add_argument('-o', help='path to labeled subset output file',
                        dest='output',
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


#============================================================================
# 2 Read files
#============================================================================
adata = sc.read_h5ad(args.input)
labels = pd.read_csv(args.labels)

print(adata.obs)
print(labels)

#============================================================================
# 3 Label
#============================================================================
# Create column that will later become index
adata.obs['later_index'] = adata.obs.index.tolist()
adata.obs['leiden'] = adata.obs['leiden'].astype(str)
labels['leiden'] = labels['leiden'].astype(str)

# Merge
adata.obs = adata.obs.merge(labels[['leiden', 'de_novo']], 
                            left_on='leiden', 
                            right_on='leiden', 
                            how='left')

# Now set 'later_index' as the index
adata.obs.set_index('later_index', inplace=True)
adata.obs.index.name = None

# Ensure that the value is a string
adata.obs['de_novo'] = adata.obs['de_novo'].astype(str)

# Check the result
print(adata)
print(adata.obs)
print(adata.obs['de_novo'].value_counts())

#============================================================================
# 4 Write file
#============================================================================
print("Writing file...")
adata.write(args.output)
