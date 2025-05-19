#!/bin/bash
#SBATCH --job-name=inference
#SBATCH --output=inference.out
#SBATCH --error=inference.err
#SBATCH --partition=GPU
#SBATCH --time=05:00:00
#SBATCH --dependency=afterok:116190

MODEL="00022-gumbel-storms-low_shot-kimg1200-color-translation-cutout"
STEP=1200
DATADIR=/soge-home/users/spet5107/mistral/alison/data/stylegan/training-runs/${MODEL}

source /lustre/soge1/users/spet5107/micromamba/etc/profile.d/micromamba.sh

micromamba activate styleGAN
python src/generate.py --outdir=${DATADIR}/results/samples --seeds=1-5000 --trunc=0.7 --network=${DATADIR}/network-snapshot-$(printf "%06d" $STEP).pkl
python src/generate_gif.py --output=${DATADIR}/results/results.gif --seed=0 --num-rows=1 --num-cols=8 --network=${DATADIR}/network-snapshot-$(printf "%06d" $STEP).pkl
python src/style_mixing.py --outdir=${DATADIR}/results/stylemixing --trunc=0.7

micromamba activate hazGAN-torch
python inference.py --model=${MODEL} --step=${STEP} 
