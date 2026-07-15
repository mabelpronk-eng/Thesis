# ============================================================
# Script: Summarize and Visualize Focal Copy Number Aberrations
# ============================================================
#
# Description:
# This script processes gene-level GISTIC copy number data for GBM samples,
# filters the cohort, categorizes copy number values into biological states
# (Deletion, Neutral, Amplification), and generates both visual and tabular outputs.
#
# Workflow:
# 1. Load GISTIC-derived gene-level copy number data from an Excel file.
# 2. Filter out specific samples (targets and normals) to define the final cohort.
# 3. Convert numeric copy number values into categorical statuses.
# 4. Compute per-gene frequencies of deletions, neutral states, and amplifications.
# 5. Create a stacked bar plot showing percentage distribution per gene.
# 6. Export a detailed Excel file containing both raw values and categorized statuses.
#
# Outputs:
# - A stacked bar plot summarizing focal aberration frequencies across genes.
# - An Excel file with per-sample gene status annotations and raw values.
#
# Purpose:
# To provide a clear overview of the prevalence of key focal genomic aberrations
# in the dataset and support downstream interpretation and reporting.
#
# Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import os

# --- 1. SETUP PATHS ---
file_path = '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/GISTIC/GISTIC/Other_focal_aberrations_GISTIC.xlsx'
save_dir = "/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/GISTIC/GISTIC/focal_aberrations/"
save_name = "focal_aberrations_summary_stacked.png"
# Define the excel output path here
output_path = os.path.join(save_dir, "Focal_Aberrations_with_Status.xlsx")

# Load the data
df = pd.read_excel(file_path)

# --- 2. CLEANING & FILTERING ---
# IDs to exclude (Targets + Normals)
targets = ['06-0221', '14-0736', '14-1402', '19-0957'] 
exclude_normals = ['06-0139', '06-0178']
all_to_exclude = targets + exclude_normals

# Filter out the samples
df_filtered = df[~df['Sample'].str.contains('|'.join(all_to_exclude))].copy()

print(f"Original sample count: {len(df)}")
print(f"Final analytical cohort size: {len(df_filtered)}")

# --- 3. CATEGORIZATION LOGIC ---
def categorize_gistic(val):
    if val < 0: return 'Deletion'
    elif val == 0: return 'Neutral'
    else: return 'Amplification'

genes_of_interest = ["MDM4", "PDGFRA", "EGFR", "PTEN", "MGMT", "CDK4", "MDM2"]

# --- 4. GENERATE STACKED BAR PLOT ---
# Prepare data for plotting
df_cat_plot = df_filtered[genes_of_interest].applymap(categorize_gistic)
counts = df_cat_plot.apply(lambda x: x.value_counts()).fillna(0).T
order = ['Deletion', 'Neutral', 'Amplification']
counts = counts.reindex(columns=order)
counts_pct = counts.div(counts.sum(axis=1), axis=0) * 100

fig, ax = plt.subplots(figsize=(12, 8))
colors = ['#3498db', '#ecf0f1', '#e74c3c'] # Blue, Grey, Red

counts_pct.plot(kind='bar', stacked=True, ax=ax, color=colors, edgecolor='black', width=0.75)

# Add percentage labels
for n, x in enumerate([*counts_pct.index.values]):
    for (proportion, y_loc) in zip(counts_pct.loc[x], counts_pct.loc[x].cumsum()):
        if proportion > 5:
            plt.text(x=n, y=(y_loc - proportion/2), s=f'{proportion:.1f}%', 
                     color="black", fontsize=10, fontweight="bold", ha="center", va="center")

plt.title(f'Frequency of Focal Genomic Aberrations (N={len(df_filtered)})', fontsize=16, fontweight='bold', pad=25)
plt.ylabel('Percentage of Samples (%)', fontsize=14)
plt.xticks(rotation=0, fontsize=12)
plt.ylim(0, 100)
plt.legend(title='CNA Status', bbox_to_anchor=(1.02, 1), loc='upper left')

# Save Plot
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.tight_layout()
plt.savefig(os.path.join(save_dir, save_name), dpi=300, bbox_inches='tight')

# --- 5. EXPORT CATEGORIZED DATA TO EXCEL ---
# Create descriptive status columns
for gene in genes_of_interest:
    df_filtered[f'{gene}_Status'] = df_filtered[gene].apply(categorize_gistic)

# Organize columns: Sample -> Status Columns -> Raw Values
status_cols = [f'{g}_Status' for g in genes_of_interest]
final_order = ['Sample'] + status_cols + genes_of_interest
df_final = df_filtered[final_order]

# Save Excel
df_final.to_excel(output_path, index=False)

print(f"Success! Plot saved to: {os.path.join(save_dir, save_name)}")
print(f"Success! Excel data saved to: {output_path}")
