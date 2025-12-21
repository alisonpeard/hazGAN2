#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --output=sm.out
#SBATCH --error=sm.err
#SBATCH --partition=Short
#SBATCH --time=12:00:00

source ~/.bashrc
source activate "/soge-home/users/spet5107/micromamba/envs/snakemake"

snakemake --profile profiles/ouce_slurm rule process_all_data --unlock
snakemake --profile profiles/ouce_slurm process_all_data --rerun-incomplete
