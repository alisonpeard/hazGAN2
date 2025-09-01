import os
import zipfile
from glob import glob


def calculate_nimgs(wildcards, years_of_samples=config["nyears"]):
    with zipfile.ZipFile(os.path.join(TRAINING_DIR, "images.zip"), 'r') as zip_ref:
        img_files = [f for f in zip_ref.namelist() if f.lower().endswith(('.png', '.jpg', '.jpeg', '.npy'))]
        nimgs = len(img_files)
    nyears = YEARN - YEAR0
    freq   = nimgs / nyears
    nsamples = int(freq * years_of_samples)
    print(f"Calculated number of samples: {nsamples} based on {nimgs} images over {nyears} years.")
    return nsamples


def get_model_path(wildcards):
    checkpoint_output = checkpoints.train_stylegan.get(**wildcards).output.outdir
    model_dirs = glob(os.path.join(checkpoint_output, "*-images-low_shot-*"))
    if not model_dirs:
        raise ValueError(f"No model directories found in {checkpoint_output}")
    latest_dir = sorted(model_dirs, key=os.path.getmtime)[-1]
    model_file = os.path.join(latest_dir, f"network-snapshot-{str(KIMG).zfill(6)}.pkl")
    if not os.path.exists(model_file):
        raise ValueError(f"Model file {model_file} does not exist")
    return model_file


rule ensure_script_executable:
    """snakemake --profile profiles/cluster ensure_script_executable"""
    output:
        touch("logs/cuda_env_ready.done")
    shell:
        """
        pwd

        chmod +x workflow/scripts/cuda_env.sh
        """


checkpoint train_stylegan:
    """>>> snakemake --profile profiles/slurm train_stylegan"""
    input:
        ready="logs/cuda_env_ready.done",
        zipfile=os.path.join(TRAINING_DIR, "images.zip")
    output:
        outdir=directory(os.path.join(GENERATED_DIR, "training-output"))
    params:
        augment="color,translation,cutout",
        kimg=KIMG
    resources:
        gpus=1
    conda:
        GPUENV
    log:
        "logs/stylegan/train.log"
    shell:
        """
        source workflow/scripts/cuda_env.sh

        mkdir -p {output.outdir}
        python workflow/src/stylegan/train.py \
            --outdir={output.outdir} \
            --data={input.zipfile} \
            --gpus={resources.gpus} \
            --DiffAugment={params.augment} \
            --kimg={params.kimg} \
            &> {log}
        """


rule generate_stylegan:
    input:
        ready="logs/cuda_env_ready.done",
        network=get_model_path
    output:
        directory(os.path.join(GENERATED_DIR, "images"))
    params:
        trunc=1.0,
        nimgs=calculate_nimgs
    conda:
        GPUENV
    log:
        "logs/stylegan/generate.log"
    shell:
        """
        source workflow/scripts/cuda_env.sh

        python workflow/src/stylegan/generate.py \
            --outdir={output} \
            --seeds=1-{params.nimgs} \
            --trunc={params.trunc} \
            --network={input.network} \
            &> {log}
        """