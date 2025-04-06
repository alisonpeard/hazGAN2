"""
Get data from SoGE cluster. Run from repository root.

snakemake --profile profiles/cluster/ all_slices --use-conda -n
snakemake --profile profiles/local/ all_slices --use-conda -n
"""
rule all_slices:
    input:
        expand(
            "{outdir}/{year}",
            outdir=PROCESSING_DIR,
            year=YEARS,
        )

rule get_data:
    output:
        netcdf=os.path.join(PROCESSING_DIR, "{year}.nc")
    params:
        year="{year}",
        indir=INDIR,
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        fields=FIELDS
    resources:
        partition="Short",
        nodes=4,
        cpus=4,
        time="00:30:00",
        job_name="era5",
        slurm="--output=sbatch_dump/era5_%A_%a.out --error=sbatch_dump/era5_%A_%a.err"
    conda:
        "environments/{env}.yaml" # option to change to path existing env
    log:
        "logs/get_data_{year}.log"
    script:
        "../scripts/get_data.py"
