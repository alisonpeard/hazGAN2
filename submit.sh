#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --output=sm.out
#SBATCH --error=sm.err
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --account=ouce-opsis

source ~/.bashrc

# Debug: Print environment info
echo "Before loading Anaconda3:"
echo "SNAKEMAKE variable: $SNAKEMAKE"
echo "PATH: $PATH"
echo ""

module load Anaconda3

echo "After loading Anaconda3:"
echo "SNAKEMAKE variable: $SNAKEMAKE"
echo "CONDA_PREFIX: $CONDA_PREFIX"
which conda
echo ""

source activate $SNAKEMAKE

echo "After activating environment:"
echo "Current conda env: $CONDA_DEFAULT_ENV"
which snakemake
echo ""

snakemake --profile profiles/slurm generate_stylegan --unlock
snakemake --profile profiles/slurm generate_stylegan --rerun-incomplete
