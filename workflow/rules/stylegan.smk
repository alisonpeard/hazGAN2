# how to make these check dir contents etc?
rule train_stylegan:
    input:
        zipfile=os.path.join(TRAINING_DIR, "jpegs.zip")
    output:
        network=os.path.join(GENERATED_DIR, "network-snapshot-000000.pkl"),
    params:
        augment="color,translation,cutout"
    conda:
        GPUENV
    log:
        file="logs/stylegan_train.log"
    shell:
        "python ../../packages/styleGAN2-DA/src/train.py --outdir={GENERATED_DIR} --data={input.zipfile} --gpus=1 --DiffAugment={params.augment} --dry-run"


rule generate_stylegan:
    input: # all .jpeg files in
        # tbc
    output:
        datadir=os.path.join(GENERATED_DIR)
    params:
        trunc=1.0
        nimgs=1000
    log:
        file="logs/stylegan_generate.log"
    shell:
        "python ../../packages/styleGAN2-DA/src/generate.py --outdir={output.datadir} --seeds=1-{params.nimgs} --trunc={params.trunc} --network=${DATADIR}/network-snapshot-$(printf "%06d" $STEP).pkl"