"""
===========================================================================
Deconvolution Performance Analysis Script
===========================================================================

This script evaluates the performance of cell type deconvolution by 
comparing true cell fractions to predicted/deconvolved fractions.

It performs the following analyses:

1. Pearson Correlation (PCC)
   - Calculates Pearson correlation for each cell type between true 
     and predicted fractions.
   - Generates a barplot of correlations per cell type.
   - Computes and prints the average correlation across all cell types.

2. Root Mean Squared Deviation (RMSD)
   - Computes RMSD per sample for each cell type.
   - Computes average RMSD per cell type and overall average RMSD.
   - Generates a boxplot showing RMSD distributions per cell type.

3. Scatterplots of True vs Predicted Fractions
   - Creates a full-scale scatterplot (0-1) for all cell types.
   - Creates a zoomed-in scatterplot (0-0.3) to highlight small fractions.
   - Identity line (y=x) is plotted for reference.

Inputs:
- CSV files containing true and deconvolved cell fractions (rows=samples, columns=cell types).

Outputs:
- Pearson correlation barplot PNG.
- RMSD boxplot PNG.
- Combined scatterplots (full-scale and zoomed-in) PNGs.

Author: Mabel Pronk
"""



#-------------------------------------------------------------------------------
# 0. Import packages and configure plotting
#-------------------------------------------------------------------------------
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error

# Global plot settings
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})

#-------------------------------------------------------------------------------
# 1. Pearson Correlation (PCC) per cell type
#-------------------------------------------------------------------------------
def correlation_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_corr):
    """Calculate Pearson correlation per cell type and generate barplot."""
    
    # Load data
    true_fractions = pd.read_csv(true_fractions_file, index_col=0).rename_axis(None)
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0)

    print("True fractions shape:", true_fractions.shape)
    print("Deconvolved fractions shape:", deconvolved_fractions.shape)
    assert true_fractions.shape == deconvolved_fractions.shape, "Shape mismatch!"

    # Dictionary to store correlations
    correlations = {}

    # Loop through each cell type (column)
    for column in true_fractions.columns:
        corr, _ = pearsonr(true_fractions[column], deconvolved_fractions[column])
        correlations[column] = corr

    # Print correlations per cell type
    for cell_type, corr in correlations.items():
        print(f"Pearson correlation for {cell_type}: {corr:.3f}")

    # Overall average correlation
    avg_corr = np.mean(list(correlations.values()))
    print(f"\nAverage Pearson Correlation across all cell types: {avg_corr:.3f}")

    # Plotting
    colors = sns.color_palette("tab10", len(correlations))
    plt.figure(figsize=(10, 6))
    plt.bar(correlations.keys(), correlations.values(), color=colors)
    plt.xlabel('Cell Types')
    plt.ylabel('Pearson Correlation')
    plt.title('Pearson Correlation between True and Deconvolved Cell Fractions')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_corr)
    print("Correlation barplot saved at:", output_corr)
    plt.close()


#-------------------------------------------------------------------------------
# 2. RMSD Boxplot per cell type
#-------------------------------------------------------------------------------
def box_RMSD_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_path):
    """Calculate RMSD per sample and per cell type, generate boxplot."""

    # Load data
    true_fractions = pd.read_csv(true_fractions_file, index_col=0).rename_axis(None)
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0)
    assert true_fractions.shape == deconvolved_fractions.shape, "Shape mismatch!"

    # DataFrame for per-sample RMSD
    rmsd_values = pd.DataFrame(index=true_fractions.index)
    rmsd_scores = {}  # Store average RMSD per cell type

    # Compute RMSD
    for cell_type in true_fractions.columns:
        rmsd_values[cell_type] = np.sqrt((true_fractions[cell_type] - deconvolved_fractions[cell_type]) ** 2)
        rmsd = np.sqrt(mean_squared_error(true_fractions[cell_type], deconvolved_fractions[cell_type]))
        rmsd_scores[cell_type] = rmsd
        print(f"Average RMSD for {cell_type}: {rmsd:.3f}")

    # Overall average RMSD
    avg_rmsd = np.mean(list(rmsd_scores.values()))
    print(f"\nOverall Average RMSD across all cell types: {avg_rmsd:.3f}")

    # Prepare data for seaborn boxplot
    rmsd_melted = rmsd_values.melt(var_name='Cell Type', value_name='RMSD')

    # Boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Cell Type', y='RMSD', data=rmsd_melted, palette='tab10', width=0.6)
    plt.title('Boxplot of RMSD for Each Cell Type')
    plt.xlabel('Cell Type')
    plt.ylabel('RMSD')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"RMSD boxplot saved at: {output_path}")
    plt.show()


#-------------------------------------------------------------------------------
# 3. Combined Scatterplot of True vs Predicted Fractions
#-------------------------------------------------------------------------------
def combined_scatter_simple(true_fractions_file, deconvolved_fractions_file,
                            output_path, output_path_zoomedin):
    """Scatterplot of true vs predicted fractions, full-scale and zoomed-in."""
    
    # Load and sort
    true_fractions = pd.read_csv(true_fractions_file, index_col=0).sort_index()
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0).sort_index()

    # Convert to long format
    df_true_long = true_fractions.reset_index().melt(
        id_vars=true_fractions.index.name or 'index', 
        var_name='cell_type', 
        value_name='true_fraction'
    )
    df_pred_long = deconvolved_fractions.reset_index().melt(
        id_vars=deconvolved_fractions.index.name or 'index',
        var_name='cell_type',
        value_name='predicted_fraction'
    )

    # Rename index column to 'sample' and merge
    df_true_long = df_true_long.rename(columns={df_true_long.columns[0]: 'sample'})
    df_pred_long = df_pred_long.rename(columns={df_pred_long.columns[0]: 'sample'})
    df = pd.merge(df_true_long, df_pred_long, on=['sample', 'cell_type'])

    # Color palette
    cell_types = df['cell_type'].unique()
    palette = sns.color_palette("tab10", len(cell_types))
    color_dict = dict(zip(cell_types, palette))

    # --- Full-scale scatterplot ---
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='true_fraction', y='predicted_fraction', 
                    hue='cell_type', palette=color_dict, alpha=0.8)
    plt.plot([0, 1], [0, 1], 'r--', label='Identity Line')
    plt.title('Scatterplot of True vs Predicted Cell Fractions')
    plt.xlabel('True Cell Fractions')
    plt.ylabel('Predicted Cell Fractions')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path)
    print("Full-scale scatterplot saved at:", output_path)
    plt.close()

    # --- Zoomed-in scatterplot (0-0.3) ---
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='true_fraction', y='predicted_fraction', 
                    hue='cell_type', palette=color_dict, alpha=0.8)
    plt.plot([0, 0.3], [0, 0.3], 'r--', label='Identity Line')
    plt.title('Zoomed-in Scatterplot (0-0.3)')
    plt.xlabel('True Cell Fractions')
    plt.ylabel('Predicted Cell Fractions')
    plt.xlim(0, 0.3)
    plt.ylim(0, 0.3)
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path_zoomedin)
    print("Zoomed-in scatterplot saved at:", output_path_zoomedin)
    plt.close()


#-------------------------------------------------------------------------------
# Main execution
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    true_fractions_file = "/net/beegfs/users/P086608/pseudobulk/level2_celltype/cell_type_fractions.csv"
    deconvolved_fractions_file = "/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/Output/fractions3.csv"

    # Pearson correlation barplot
    output_corr = "/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/Output/correlation_barplot.png"
    correlation_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_corr)

    # RMSD boxplot
    output_rmsd = "/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/Output/rmsd_boxplot.png"
    box_RMSD_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_rmsd)

    # Combined scatterplots
    output_scatter = "/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/Output/combined_scatter.png"
    output_scatter_zoom = "/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/Output/combined_scatter_zoom.png"
    combined_scatter_simple(true_fractions_file, deconvolved_fractions_file, output_scatter, output_scatter_zoom)
