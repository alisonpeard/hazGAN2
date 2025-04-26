checkpoint train_stylegan:
    input:
        zipfile=os.path.join(TRAINING_DIR, "images.zip")
    output:
        directory(os.path.join(GENERATED_DIR, "training-output"))
    params:
        augment="color,translation,cutout",
        kimg=KIMG
    conda:
        GPUENV
    log:
        "logs/stylegan_train.log"
    shell:
        """
        mkdir -p {params.outdir}
        python ../../packages/styleGAN2-DA/src/train.py \
            --outdir={params.outdir} \
            --data={input.zipfile} \
            --gpus=1 \
            --DiffAugment={params.augment} \
            --kimg={params.kimg} \
            &> {log}
        """


rule generate_stylegan:
    input:
        network=os.path.join(TRAINING_DIR,
                            "training-output",
                            "00001-images-low_shot-color-translation-cutout",
                            "network-snapshot-{KIMG}.pkl")
    output:
        directory(os.path.join(GENERATED_DIR, "images"))
    params:
        trunc=1.0
        nimgs=1000
    log:
        "logs/stylegan_generate.log"
    shell:
        """
        python ../../packages/styleGAN2-DA/src/generate.py \
            --outdir={output} \
            --seeds=1-{params.nimgs} \
            --trunc={params.trunc} \
            --network={input.network} \
            &> {log}
        """