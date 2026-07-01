# ============================================================
# Script: Absolute Copy Number Profile Visualization
# ============================================================
#
# Description:
# This script generates per-sample absolute copy number profiles
# using segmented copy number data and ACE model fits (cellularity + ploidy).
# It visualizes tumor copy number architecture corrected for purity and ploidy.
#
# Workflow:
# 1. Load ACE helper functions
# 2. Iterate over segmented copy number files
# 3. Convert segments to copy number templates
# 4. Retrieve best-fit purity/ploidy estimates per sample
# 5. Generate absolute CN profiles using ACE
# 6. Save plots per sample
#
# Output:
# - PNG plots of absolute copy number profiles per tumor sample
#
# ============================================================

# Load required functions
source('/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/scripts/ACE_functions.R')


# Paths
segments_dir <- '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/copynumber/segments/'
fits_dir <- '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/copynumber/ACE/'
output_dir <- '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/plots/Copy_number_profile_absolute/'

# Create output directory if it doesn't exist
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# List all segment files
# Each file corresponds to one tumor sample segmentation output
segment_files <- list.files(segments_dir, pattern = "_segments.txt$", full.names = TRUE)

# Loop over files
for(seg_file in segment_files) {
  
  # Read segments
  Segments <- read.delim(seg_file, stringsAsFactors = FALSE)
  
  # Generate template
  # Converts segmented CN data into genome-wide template format (log2 scale)
  template <- segmentstotemplate(Segments, log = 2)
  
  # Determine sample name from file name
  sample_name <- gsub("_segments.txt", "", basename(seg_file))
  
  # Read bestfit ACE model results for this sample
  bestfit_file <- file.path(fits_dir, sample_name, "squaremodel", "fits.txt")

  # Skip sample if ACE fit is missing
  if(!file.exists(bestfit_file)) {
    warning(paste("Bestfit file not found for", sample_name, "- skipping."))
    next
  }

  bestfit <- read.delim(bestfit_file, stringsAsFactors = FALSE)
  
  # Use first row (best model fit)
  # Contains estimated tumor purity (cellularity) and ploidy
  best_row <- bestfit[1, ]
  cellularity <- best_row$cellularity
  ploidy <- best_row$ploidy
  
  # Generate absolute copy number plot 
  #NOT ADDED -> error and standard could add these to 
  p <- singleplot(template, cellularity = cellularity, ploidy = ploidy,
                  title = paste0(sample_name, " CN profile"))
  
  # Save plot
  output_file <- file.path(output_dir, paste0(sample_name, "_CN_profile.png"))
  ggsave(output_file, plot = p, width = 8, height = 4)
  
  # Log progress
  cat("Saved plot for", sample_name, "\n")
}
