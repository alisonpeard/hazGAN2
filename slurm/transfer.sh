# edit this as needed
INDIR="/soge-home/users/spet5107/code/hazGAN2/projects/poweruk_winter/results/training"
INLOC="linux.ouce.ox.ac.uk"
INUSER="spet5107"

OUTDIR="/Users/alison/Documents/dphil/data/hazGAN2/projects/poweruk_winter/results/training"

scp -r ${INUSER}@${INLOC}:${INDIR} ${OUTDIR}