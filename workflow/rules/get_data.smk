"""
Get data from SoGE cluster. Run from repository root.

```bash
micromamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge conda=24.7.1
conda config --set channel_priority strict
snakemake --profile profiles/cluster/ all_slices --use-conda -n 
```
"""
rule all_slices:
    input:
        expand(
            os.path.join(PROCESSING_DIR, "{year}.nc"),
            year=YEARS,
        )

rule get_data:
    output:
        netcdf=os.path.join(PROCESSING_DIR, "{year}.nc")
    params:
        year="{year}",
        conda=CONDA,
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
        slurm="--output=sbatch_dump/%A_%a.out --error=sbatch_dump/%A_%a.err"
    conda:
        os.path.join("..", "..", CONDA) # or existing named env since 6.14.0 (discouraged)
    log:
        "logs/get_data_{year}.log"
    script:
        "../scripts/get_data.py"
