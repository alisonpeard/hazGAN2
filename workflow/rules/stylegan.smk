import zipfile

def calculate_nimgs(wildcards, years_of_samples=1000):
    with zipfile.ZipFile(os.path.join(TRAINING_DIR, "images.zip"), 'r') as zip_ref:
        img_files = [f for f in zip_ref.namelist() if f.lower().endswith(('.png', '.jpg', '.jpeg'))]
        nimgs = len(img_files)
    nyears = YEARN - YEAR0
    freq   = nimgs / nyears
    nsamples = int(nimgs * years_of_samples)
    return nsamples
    

checkpoint train_stylegan:
    input:
        zipfile=os.path.join(TRAINING_DIR, "images.zip")
    output:
        outdir=directory(os.path.join(GENERATED_DIR, "training-output"))
    params:
        augment="color,translation,cutout",
        kimg=KIMG
    conda:
        GPUENV
    log:
        "logs/stylegan/train.log"
    shell:
        """
        mkdir -p {output.outdir}
        python ../../packages/styleGAN2-DA/src/train.py \
            --outdir={output.outdir} \
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
                            f"network-snapshot-{KIMG}.pkl")
    output:
        directory(os.path.join(GENERATED_DIR, "images"))
    params:
        trunc=1.0,
        nimgs=1000 #! want to automate this based on number of images in training set (images.zip) and YEARN - YEAR0
    conda:
        GPUENV
    log:
        "logs/stylegan/generate.log"
    shell:
        """
        python ../../packages/styleGAN2-DA/src/generate.py \
            --outdir={output} \
            --seeds=1-{params.nimgs} \
            --trunc={params.trunc} \
            --network={input.network} \
            &> {log}
        """