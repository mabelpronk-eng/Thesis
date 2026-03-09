#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# compile_transcriptome_matrix.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# Compile transcriptome matrices to one big matrix
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl) & Mabel Pronk (m.pronk3@amsterdamumc.nl)
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(dplyr)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.0 Read data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("Reading data...")
project = "TCGA-LGG"
# read table
RNA_sampletable <- read.delim("/net/beegfs/users/P086608/bulkRNA_glioma/data/LGG/sampletable/RNA_Sampletable.tsv", stringsAsFactors = F) %>%
    # obtain tumors
    filter(sample_type == "Primary Tumor") %>%
    # fetch correct ids
    mutate(TCGA_sample = purrr::map_chr(cases,~paste(strsplit(.x,"-")[[1]][1:3],collapse = "-")),
           project_short = gsub("TCGA-","",project),
           filepath = paste0("/net/beegfs/users/P086608/bulkRNA_glioma/data/LGG/TCGA-LGG/TCGA-LGG/Transcriptome_Profiling/Gene_Expression_Quantification/",id,"/",file_name))

print(paste0("Found ",nrow(RNA_sampletable)," samples."))
Read_expression_data <- function(file){
    counts <- data.table::fread(file, header = TRUE) %>%
        filter(startsWith(gene_id, "ENSG")) %>%
        select(gene = gene_id, gene_name = gene_name, count = unstranded)
    return(counts)
}

# Read all expression data
print("Reading expression data...")
transcriptome_data <-
    tibble::tribble(
            ~project,~sample,~data,
            RNA_sampletable$project_short,RNA_sampletable$TCGA_sample,purrr::map(RNA_sampletable$filepath,Read_expression_data)) %>%
    tidyr::unnest()%>%
    tidyr::unnest()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.0 Reformat data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("Reformatting data...")
transcriptome_data_disease <-
    transcriptome_data %>%
    filter(project == gsub("TCGA-","",project)) %>%
    select(-project) %>%
    unique()

# Retrieve wide format
print("Converting to wide format...")
transcriptome_data_disease <-
    transcriptome_data_disease %>%
    tidyr::pivot_wider(id_cols = c(gene, gene_name), names_from = sample, values_from = count, values_fn = sum)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3.0 write to file
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("Writing to file...")
write.table(transcriptome_data_disease, paste0("data/LGG/",project,"/Transcriptome_matrix.txt"), quote = F, row.names = F,sep = "\t")
write.table(transcriptome_data_disease, paste0("data/LGG/",project,"/Transcriptome_matrix.csv"), quote = F, row.names = F,sep = ",")
