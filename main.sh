#!/bin/bash
#SBATCH --job-name=manager
#SBATCH --nodes=4
#SBATCH --cpus-per-task=4
#SBATCH --output=sbatch_dump/bigsnake_%A.out
#SBATCH --error=sbatch_dump/bigsnake_%A.err
#SBATCH --partition=Medium

# source activate micromamba
# micromamba activate snakemake

snakemake --profile profiles/slurm/ process_all_data

# activate snakemake
# sbatch manager.sh