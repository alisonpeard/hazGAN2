#!/bin/bash
#SBATCH --job-name=main
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=sbatch_dump/main_%A.out
#SBATCH --error=sbatch_dump/main_%A.err
#SBATCH --partition=Medium

# might need to run these before:
# snakemake --profile profiles/slurm/ process_all_data --unlock
# bash cleanrepo.sh

snakemake --profile profiles/slurm/ process_all_data