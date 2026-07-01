#!/bin/bash
#SBATCH --job-name=DEMO_statescope_cpu_pseudobulk
#SBATCH --partition=defq                 # adjust for your cluster 
#SBATCH --cpus-per-task=32              # total CPU cores
#SBATCH --mem=500G  
#SBATCH --time=12:00:00
#SBATCH --output=/net/beegfs/users/P086608/Deconvolution_pipeline-not-used-yet/Output/logs/%x_%j.out
#SBATCH --error=/net/beegfs/users/P086608/Deconvolution_pipeline-not-used-yet/Output/logs/%x_%j.err


# Activate environment 
source /home/P086608/miniconda3/etc/profile.d/conda.sh 

# activate your env 
conda activate /home/P086608/miniconda3/envs/snakemake 

cd /net/beegfs/users/P086608/Deconvolution_pipeline-not-used-yet

# Run (CPU)
snakemake -s Snakefile --use-conda -c 40
