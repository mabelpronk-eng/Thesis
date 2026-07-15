################################################################################
# Generate copy number profile
################################################################################
# This script generates a multi-panel copy number visualization for a GBM sample.
#
# It combines:
#   Panel A → Genome-wide copy number profile across all autosomes (1–22)
#   Panel B → Chromosome 9 overview (showing CN segments + key genes)
#   Panel C → Zoomed-in view of chromosome 9p21.3 region (CDKN2A/B, IFNE locus)
#
# Input:
#   - Segmentation file from copy number analysis (per sample)
#   - Contains chromosome segments and log2 copy number ratios
#
# Output:
#   - Combined PNG figure with three panels (A–C)
#
# Key biological focus:
#   - Chromosome 9 alterations in GBM (including 9p21.3 deletion region)
#   - Genes of interest: IFNE, CDKN2A, CDKN2B
#
# Dependencies:
#   - patchwork (plot assembly)
#   - ggplot2 (visualisation)
#   - dplyr (data wrangling)
#   - ggrepel (non-overlapping labels)
#
# Author: Mabel Pronk (m.pronk3@amsterdamumc.nl)
################################################################################


library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)

# --- 1. SETTINGS & PATHS ---
# Define sample ID to analyse
sample_id <- "TCGA-06-0645" 

# Path to segmented copy number data for this sample
seg_file  <- paste0("/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/copynumber/segments/", sample_id, "_segments.txt")

# Output file name for final combined plot
out_file  <- paste0(sample_id, "_fixed_combined_plot.png")

# Load helper functions for segment processing (external script)
source('/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/scripts/ACE_functions.R')

# --- 2. DEFINE FUNCTIONS ---

# ------------------------------------------------------------------------------
# Panel A: Genome-wide copy number plot
# ------------------------------------------------------------------------------
# Creates a whole-genome view of copy number values across chromosomes 1–22.
# Each point represents a segment log2 ratio value positioned along the genome.
plot_template <- function(template, title = "Copy Number Profile", subtitle = NULL) {
  
  df <- template
  
  # Convert chromosome X to numeric (23) for consistent ordering
  df$chr <- as.numeric(gsub("X", "23", as.character(df$chr)))
  
  # Keep only autosomes (1–22)
  df <- df %>% filter(chr %in% 1:22)
  
  # Compute chromosome sizes for plotting cumulative genome position
  chr_sizes <- df %>% 
    group_by(chr) %>% 
    summarise(chr_len = max(end)) %>% 
    arrange(chr)
  
  chr_sizes$cumstart <- c(0, cumsum(head(chr_sizes$chr_len, -1)))
  
  # Map segment positions into a continuous genome coordinate
  df <- df %>% 
    left_join(chr_sizes, by = "chr") %>% 
    mutate(genome_pos = start + cumstart)
  
  # Plot genome-wide CN profile
  ggplot(df) +
    geom_point(aes(x = genome_pos, y = copynumbers), size = 0.3, alpha = 0.6) +
    
    # Add chromosome boundaries
    geom_vline(xintercept = chr_sizes$cumstart, color = "grey70") +
    
    # Chromosome labels centered
    scale_x_continuous(
      breaks = chr_sizes$cumstart + chr_sizes$chr_len/2, 
      labels = as.character(1:22),
      expand = c(0,0)
    ) +
    
    theme_bw() +
    
    # Fixed y-scale for comparability
    scale_y_continuous(breaks = seq(-6, 6, 2), limits = c(-6, 6)) +
    
    labs(
      title = title,
      subtitle = subtitle,
      x = "Chromosomes",
      y = expression(log[2]~ratio)
    )
}

# ------------------------------------------------------------------------------
# Panels B & C: Chromosome 9 visualization
# ------------------------------------------------------------------------------
# Creates CN segment plot for chromosome 9 with optional zoom window.
# Highlights key GBM genes: IFNE, CDKN2A, CDKN2B.
make_chr9_plot <- function(df, xlim = NULL, subtitle = NULL, highlight_zoom = NULL, is_zoom = FALSE) {
  
  # Define genes of interest on chromosome 9 (Mb positions)
  genes <- data.frame(
    name = c("IFNE", "CDKN2A", "CDKN2B"),
    pos = c(21.481, 21.968, 22.003) 
  )
  
  # Base plot
  p <- ggplot(df) +
    theme_classic(base_size = 11) +
    geom_hline(yintercept = 0, color = "gray60")
  
  # Highlight region of interest (e.g. 9p21.3)
  if(!is.null(highlight_zoom)) {
    p <- p + 
      annotate("rect", xmin = highlight_zoom[1], xmax = highlight_zoom[2], 
               ymin = -Inf, ymax = Inf, fill = "royalblue", alpha = 0.1) +
      geom_vline(xintercept = highlight_zoom, color = "royalblue", linetype = "dotted")
  }
  
  # Optional background shading for zoom panel
  if(is_zoom) {
    p <- p + theme(
      panel.background = element_rect(fill = alpha("royalblue", 0.1), color = NA)
    )
  }

  # Add CN segment lines
  p <- p +
    geom_segment(
      aes(x = Start_Mbp, xend = End_Mbp, y = Segment_Mean, yend = Segment_Mean),
      color = 'black', linewidth = 1.2, lineend = "round"
    ) +
    
    # Add gene positions
    geom_vline(
      data = genes,
      aes(xintercept = pos),
      color = "darkblue",
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.8
    ) +
    
    # Gene labels
    geom_text_repel(
      data = genes,
      aes(x = pos, y = 1.6, label = name),
      fontface = "italic",
      size = 3.5,
      nudge_y = 0.2
    ) +
    
    labs(
      subtitle = subtitle,
      x = "Position on Chromosome 9 (Mbp)",
      y = expression(Log[2]~Ratio)
    ) +
    
    coord_cartesian(ylim = c(-4, 2), xlim = xlim, expand = FALSE) +
    
    theme(
      legend.position = "none",
      panel.border = element_rect(fill=NA, colour = "black")
    )
  
  return(p)
}

# --- 3. RUN PIPELINE ---
# Only run if segmentation file exists
if (file.exists(seg_file)) {
  
  # Read copy number segments
  Segments <- read.delim(seg_file, stringsAsFactors = FALSE)
  
  # Convert segments into genome-wide plotting format (external function)
  template_data <- segmentstotemplate(Segments)
  
  # Filter chromosome 9 and convert positions to Mbp
  df_chr9 <- Segments %>%
    filter(Chromosome %in% c(9, "9", "chr9")) %>%
    mutate(Start_Mbp = Start / 1e6, End_Mbp = End / 1e6)
  
  # Create plots
  p_A <- plot_template(template_data, title = paste("Sample:", sample_id), subtitle = 'Genome-wide overview')
  p_B <- make_chr9_plot(df_chr9, subtitle = "Chromosome 9", highlight_zoom = c(15, 30))
  p_C <- make_chr9_plot(df_chr9, subtitle = "9p21.3 Zoom (15-30 Mbp)", xlim = c(15, 30), is_zoom = TRUE)
  
  # Combine panels into single figure
  final_plot <- p_A | p_B | p_C + 
    plot_annotation(tag_levels = 'A')
  
  # Save final figure
  ggsave(
    '/net/beegfs/users/P086608/bulkRNA_glioma/data/GBM/plots/test_final/final_06-0645.png',
    plot = final_plot,
    width = 16,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  message("Success! Plot saved as: ", out_file)
}
