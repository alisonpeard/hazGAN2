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

# qsub -N PIPE -cwd -j yes python snakemake --cluster "ssh spet5107@headnode_address 'qsub -N pipe_task -j yes -cwd -S /bin/sh ' " -j

snakemake --profile profiles/slurm
