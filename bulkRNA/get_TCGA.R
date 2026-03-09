#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get_TCGA.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Fetch data from GDC data portal 
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
# Adapted by: Mabel Pronk (m.pronk3@amsterdamumc.nl)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(TCGAbiolinks)

#-------------------------------------------------------------------------------
# 0. Preparations
#-------------------------------------------------------------------------------
# Define projects 
Projects <- c("TCGA-GBM", 'TCGA-LGG')
# Determine whether to use legacy database or not (original TCGA data - older format of the data, not harmonized)
legacy <- FALSE

# Create directories
target_data_dir <- "data/" ; dir.create(target_data_dir, showWarnings = FALSE)
target_data_dir_sampletables <- "data/sampletable/" ; dir.create(target_data_dir_sampletables, showWarnings = FALSE)
target_data_dir_puritites <- "data/purity/" ; dir.create(target_data_dir_puritites, showWarnings = FALSE)
target_data_dir_clinical <- "data/clinical/" ; dir.create(target_data_dir_clinical, showWarnings = FALSE)
target_data_dir_GBM <- "data/GBM/" ; dir.create(target_data_dir_GBM , showWarnings = FALSE)
target_data_dir_LGG <- "data/LGG/" ; dir.create(target_data_dir_LGG, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 1.1  DATASET1: Create query and summary files
#-------------------------------------------------------------------------------

# Create GBM queries
query_rna_GBM <- GDCquery(project = Projects[1],
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts")
query_cna_GBM <- GDCquery(
         project = Projects[1],
         data.category = "Copy Number Variation",
         data.type = "Masked Copy Number Segment")

# Get all patients that have DNA+RNA data available
common.patients_GBM <- Reduce(intersect,
    list(substr(getResults(query_rna_GBM, cols = "cases"), 1, 12),
    substr(getResults(query_cna_GBM, cols = "cases"), 1, 12)))

# Obtain queries for overlapping patients
query_rna_GBM <- GDCquery(project = Projects[1],
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts",
         barcode = common.patients_GBM)
query_cna_GBM <- GDCquery(
         project = Projects[1],
         data.category = "Copy Number Variation",
         data.type = "Masked Copy Number Segment",
         barcode = common.patients_GBM)


results_query_rna_GBM <- getResults(query_rna_GBM)
results_query_cna_GBM <- getResults(query_cna_GBM)

# Create LGG queries
query_rna_LGG <- GDCquery(project = Projects[2],
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts")
query_cna_LGG <- GDCquery(
         project = Projects[2],
         data.category = "Copy Number Variation",
         data.type = "Masked Copy Number Segment")

# Get all patients that have DNA+RNA data available
common.patients_LGG <- Reduce(intersect,
    list(substr(getResults(query_rna_LGG, cols = "cases"), 1, 12),
    substr(getResults(query_cna_LGG, cols = "cases"), 1, 12)))

# Obtain queries for overlapping patients
query_rna_LGG <- GDCquery(project = Projects[2],
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "STAR - Counts",
         barcode = common.patients_LGG)
query_cna_LGG <- GDCquery(
         project = Projects[2],
         data.category = "Copy Number Variation",
         data.type = "Masked Copy Number Segment",
         barcode = common.patients_LGG)


results_query_rna_LGG <- getResults(query_rna_LGG)
results_query_cna_LGG <- getResults(query_cna_LGG)
#-------------------------------------------------------------------------------
# 1.3  Create sample tables
#-------------------------------------------------------------------------------

# write CNA sampletable to file
write.table(rbind(results_query_rna_GBM,results_query_rna_LGG),
            paste0(target_data_dir_sampletables,"RNA_Sampletable.tsv"),
            quote = F, row.names = F, sep = "\t")

# write CNA sampletable to file
write.table(rbind(results_query_cna_GBM,results_query_cna_LGG),
            paste0(target_data_dir_sampletables,"CNA_Sampletable.tsv"),
            quote = F, row.names = F, sep = "\t")

#-------------------------------------------------------------------------------
# Download TCGA files
GDCdownload(query_rna_GBM, method = "api", directory = target_data_dir_GBM)
print(paste0('--------------------------------Done RNA, Project: ', Projects[1], '--------------------------------'))
GDCdownload(query_cna_GBM, method = "api", directory = target_data_dir_GBM)
print(paste0('--------------------------------Done CNA, Project: ', Projects[1], '--------------------------------'))

GDCdownload(query_rna_LGG, method = "api", directory = target_data_dir_LGG)
print(paste0('--------------------------------Done RNA, Project: ', Projects[2], '--------------------------------'))
GDCdownload(query_cna_LGG, method = "api", directory = target_data_dir_LGG)
print(paste0('--------------------------------Done CNA, Project: ', Projects[2], '--------------------------------'))

# 1: Obtain tumor purities
#-------------------------------------------------------------------------------
Tumor_purities <- Tumor.purity[Tumor.purity$Cancer.type %in% c("GBM","LGG"),]
# write to file

# write to file
write.table(Tumor_purities,
            paste0(target_data_dir_puritites, "_Tumor_purities.tsv"),
            quote = F, row.names = F, sep = "\t")

#-------------------------------------------------------------------------------
# 2.2 Fetch clinical data 
#-------------------------------------------------------------------------------
Clinical_LGG <- GDCquery_clinic(Projects[1])
Clinical_GBM <- GDCquery_clinic(Projects[2])

write.table(rbind(Clinical_LGG[,colnames(Clinical_LGG) %in% colnames(Clinical_GBM)],Clinical_GBM),
            paste0(target_data_dir_clinical,"Clinical_data.tsv"),
            quote = F, row.names = F, sep = "\t")


print("Done!")

#-------------------------------------------------------------------------------
# If copy number data is too big and the download is not working, use this in an 
# interactive R session and manually adjust the number of files downloaded at once.
# Usually download works if the size of the files that are downloaded is below 4MB
# for CNA
#-------------------------------------------------------------------------------
# query_subset = query_cna
# query_subset$results[[1]]=query_cna$results[[1]][1:300,]
# GDCdownload(query_subset, method = "api", directory = target_data_dir_project)
# query_subset$results[[1]]=query_cna$results[[1]][1:600,]
# GDCdownload(query_subset, method = "api", directory = target_data_dir_project)
# query_subset$results[[1]]=query_cna$results[[1]][1:1000,]
# GDCdownload(query_subset, method = "api", directory = target_data_dir_project)
