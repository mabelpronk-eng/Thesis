# GBM Thesis Project

This repository contains the analysis scripts used for my glioblastoma thesis project.
The project involved the analysis of multiple data types, including DNA copy number data, single-cell RNA sequencing (scRNA-seq) data, and bulk RNA-seq data. In addition, a major component of the project focused on cell-type deconvolution analysis.

The scripts are organized into the following folders:

## Folders

- **CNA**
  Scripts for processing and analyzing DNA copy number data, including ACE-based tumor purity estimation, retrieval of GISTIC copy number status, and generation of copy number profiles. 

- **Single-cell RNA-seq data**  
  Included in this folder are script concerning: downloading the data, single cell atlas construction, visualization (dotplots & UMAP) of the data and pseudobulk sample generation.

- **Bulk RNA-seq data**  
  Scripts for processing and analyzing bulk transcriptomic data, including downloading TCGA datasets and preparing bulk RNA-seq data for downstream analyses.

- **Deconvolution**  
  Scripts for estimating cell-type proportions from bulk RNA-seq data using deconvolution methods (CIBERSORTx and Statescope). Includes scripts for performance evaluation, correlation analyses, visualization, and statistical comparisons.
