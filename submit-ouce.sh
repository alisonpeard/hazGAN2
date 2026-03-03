#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --output=sm.out
#SBATCH --error=sm.err
#SBATCH --partition=Short
#SBATCH --time=12:00:00

source ~/.bashrc
source activate "/soge-home/users/spet5107/micromamba/envs/snakemake"

snakemake --profile profiles/slurm get_all_years -n
snakemake --profile profiles/slurm get_all_years --unlock
snakemake --profile profiles/slurm get_all_years --rerun-incomplete
