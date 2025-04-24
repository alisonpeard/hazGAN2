rule train_stylegan:
    # input:
    #     indir=os.path.join(TRAINING_DIR, "jpegs")
    # output:
    #     netcdf=os.path.join(GENERATED_DIR, ".nc")
    params:
        year="{year}",
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        timecol=TIMECOL,
        fields=FIELDS
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/get_%A_%a.out --error=sbatch_dump/get_%A_%a.err"
    conda:
        GPUENV
    log:
        file="logs/get_{year}.log"
    script:
        "../scripts/get_data.py"