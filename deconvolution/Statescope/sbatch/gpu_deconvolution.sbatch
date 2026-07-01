#!/bin/bash
#SBATCH --job-name=statescope_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:2g.24gb
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --chdir=/net/beegfs/users/P086608/Deconvolution_pipeline_test
#SBATCH --output=/net/beegfs/users/P086608/Deconvolution_pipeline_test/Output/logs/%x_%j.out
#SBATCH --error=/net/beegfs/users/P086608/Deconvolution_pipeline_test/Output/logs/%x_%j.err
#SBATCH --mail-user=m.pronk3@amsterdamumc.nl
#SBATCH --mail-type=END,FAIL

# Activate environment
source /home/P086608/miniconda3/etc/profile.d/conda.sh

# activate your env
conda activate /home/P086608/miniconda3/envs/statescope

cd /net/beegfs/users/P086608/Deconvolution_pipeline_test

# Run (GPU)

snakemake -s Snakefile -j 1 --resources gpu=1 
