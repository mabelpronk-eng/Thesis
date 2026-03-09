#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get_TCGA.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Fetch data from GDC data portal 
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl) & Mabel Pronk (m.pronk3@amsterdamumc.nl)
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(TCGAbiolinks)

#-------------------------------------------------------------------------------
# 0. Preparations
#-------------------------------------------------------------------------------
# Define projects 

Project <- "TCGA-LGG"
# Determine whether to use legacy database or not (original TCGA data - older format of the data, not harmonized)
legacy <- FALSE

# Create directories
target_data_dir <- "data/" ; dir.create(target_data_dir, showWarnings = FALSE)
target_data_dir_sampletables <- "data/sampletable/" ; dir.create(target_data_dir_sampletables, showWarnings = FALSE)
target_data_dir_puritites <- "data/purity/" ; dir.create(target_data_dir_puritites, showWarnings = FALSE)
target_data_dir_clinical <- "data/clinical/" ; dir.create(target_data_dir_clinical, showWarnings = FALSE)
target_data_dir_project <- paste0("data/", Project, "/") ; dir.create(target_data_dir_project, showWarnings = FALSE)


#-------------------------------------------------------------------------------
# 1.1  DATASET1: Create query and summary files
#-------------------------------------------------------------------------------

# Create queries
query_rna <- GDCquery(project = Project,
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts")
query_cna <- GDCquery(
         project = Project,
         data.category = "Copy Number Variation",
         data.type = "Masked Copy Number Segment")

# Get all patients that have DNA+RNA data available
common.patients <- Reduce(intersect,
    list(substr(getResults(query_rna, cols = "cases"), 1, 12),
    substr(getResults(query_cna, cols = "cases"), 1, 12)))

# Obtain queries for overlapping patients
query_rna <- GDCquery(project = Project,
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts",
         barcode = common.patients)
query_cna <- GDCquery(
         project = Project,
         data.category = "Copy Number Variation",
         data.type = "Masked Copy Number Segment",
         barcode = common.patients)


results_query_rna <- getResults(query_rna)
results_query_cna <- getResults(query_cna)


#-------------------------------------------------------------------------------
# 1.3  Create sample tables
#-------------------------------------------------------------------------------
# write RNA sampletable to file
write.table(results_query_rna,
            paste0(target_data_dir_sampletables,"RNA_Sampletable.tsv"),
            quote = F, row.names = F, sep = "\t")

# write CNA sampletable to file
write.table(results_query_cna,
            paste0(target_data_dir_sampletables,"CNA_Sampletable.tsv"),
            quote = F, row.names = F, sep = "\t")

#-------------------------------------------------------------------------------
# Download TCGA files
GDCdownload(query_rna, method = "api", directory = target_data_dir_project)
print(paste0('--------------------------------Done RNA, Project: ', Project, '--------------------------------'))
GDCdownload(query_cna, method = "api", directory = target_data_dir_project, files.per.chunk = 100)
print(paste0('--------------------------------Done CNA, Project: ', Project, '--------------------------------'))

# 1: Obtain tumor purities
#-------------------------------------------------------------------------------
disease <- strsplit(Project, "-")[[1]][2]
Tumor_purities <- Tumor.purity[Tumor.purity$Cancer.type %in% disease,]
# write to file
write.table(Tumor_purities,
            paste0(target_data_dir_puritites, disease, "_Tumor_purities.tsv"),
            quote = F, row.names = F, sep = "\t")

#-------------------------------------------------------------------------------
# 2.2 Fetch clinical data 
#-------------------------------------------------------------------------------

Clinical <- GDCquery_clinic(Project)

library(dplyr)

# Flatten any list-columns
Clinical_flat <- Clinical %>%
  mutate(across(where(is.list), ~ sapply(., paste, collapse = ";")))

write.table(Clinical_flat,
            paste0(target_data_dir_clinical, disease ,"_Clinical_data.tsv"),
            quote = F, row.names = F, sep = "\t")

print("Done!")

