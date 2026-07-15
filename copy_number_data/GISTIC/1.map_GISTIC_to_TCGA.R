################################################################################
# Obtain GISTIC values (map gistic to TCGA samples)
################################################################################
# This script integrates GISTIC copy number alteration data with local GBM
# segmentation sample information to extract and summarise focal genomic
# aberrations in selected cancer-related genes.
#
# Workflow overview:
#
# 1. Reads available CNV segment sample IDs from local segmentation files
# 2. Loads gene-level GISTIC2 copy number thresholded data
# 3. Harmonises sample IDs between GISTIC output and local segment dataset
# 4. Filters for genes of interest (e.g. EGFR, PTEN, CDK4, MDM2, etc.)
# 5. Retains only samples present in the segmentation dataset
# 6. Transposes the matrix so samples become rows (for downstream analysis)
# 7. Exports the final curated gene-by-sample CNV matrix to Excel
#
# Output:
#   - Excel file containing selected gene-level copy number alterations per sample
#
# Purpose:
#   - To enable comparison of focal CNV events across key GBM driver genes
#   - To link GISTIC calls with sample-level CNV segmentation data
#
# Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
################################################################################

library(dplyr)
library(readr)
library(openxlsx)

#------------------------------------------------------------
# 1. Paths
#------------------------------------------------------------
segment_dir <- "/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/copynumber/segments"
gistic_file <- "/net/beegfs/users/P086608/gdac.broadinstitute.org_GBM-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt"
out_file <- "/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/classification/GISTIC/Other_focal_aberrations_GISTIC.xlsx"

#------------------------------------------------------------
# 2. Get sample names from segment files
#------------------------------------------------------------
segment_files <- list.files(segment_dir, pattern = "_segments.txt$", full.names = FALSE)

samples <- gsub("_segments.txt$", "", segment_files)

length(samples)   # sanity check

#------------------------------------------------------------
# 3. Read GISTIC gene-level file
#------------------------------------------------------------
GBM <- read.delim(
  gistic_file,
  check.names = FALSE,
  stringsAsFactors = FALSE)

#------------------------------------------------------------
# 3b. Harmonize GISTIC sample IDs to match segment samples
#------------------------------------------------------------
gistic_samples <- colnames(GBM)[-(1:4)]  # skip metadata columns

# Truncate to TCGA-XX-XXXX-01A format
gistic_samples_short <- substr(gistic_samples, 1, 12)

colnames(GBM)[-(1:4)] <- gistic_samples_short

#------------------------------------------------------------
# 4. Select genes of interest
#------------------------------------------------------------
genes_of_interest <- c("MDM4", 'PTEN', 'MGMT', 'MDM2','CDK4', "PDGFRA", "EGFR")

GBM_sel <- GBM %>%
  filter(`Gene Symbol` %in% genes_of_interest)

#------------------------------------------------------------
# 5. Keep only samples that exist in segments folder
#------------------------------------------------------------
sample_cols <- intersect(samples, colnames(GBM_sel))

GBM_sel <- GBM_sel %>%
  select(`Gene Symbol`, all_of(sample_cols))

#------------------------------------------------------------
# 6. Transpose: samples as rows
#------------------------------------------------------------
GBM_mat <- as.data.frame(t(GBM_sel[,-1]))
colnames(GBM_mat) <- GBM_sel$`Gene Symbol`

GBM_mat$Sample <- rownames(GBM_mat)

GBM_final <- GBM_mat %>%
  relocate(Sample)

#------------------------------------------------------------
# 7. Write to Excel
#------------------------------------------------------------
# Create the directory if it doesn't exist
if (!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}

# Write the file
write.xlsx(
  GBM_final,
  file = out_file,
  overwrite = TRUE
)

cat("Excel file successfully written to:\n", out_file, "\n")
