# how to make these check dir contents etc?
rule train_stylegan:
    input:
        os.path.join(TRAINING_DIR, "jpegs.zip")
    output:
        GENERATED_DIR
    params:
        year="{year}",
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        timecol=TIMECOL,
        fields=FIELDS
    conda:
        GPUENV
    log:
        file="logs/stylegan_{year}.log"
    shell:
        "python ../../packages/styleGAN2-DA/src/train.py --outdir={GENERATED_DIR} --data={input} --gpus=1 --DiffAugment=color,translation,cutout -dry-run"


rule generate_stylegan: