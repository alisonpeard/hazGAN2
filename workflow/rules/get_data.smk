"""
Get data from SoGE cluster. Run from repository root.

```bash
# environment set up
micromamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge conda=24.7.1
install snakemake-executor-plugin-slurm

# run rule project root and login node
micromamba activate snakemake
snakemake --profile profiles/cluster/ --executor slurm get_data --use-conda
```
"""
rule get_data:
    """Rule to process all years for the project."""
    input:
        expand(
            os.path.join(PROCESSING_DIR, "daily_{year}.nc"),
            year=YEARS
        )

rule get_year:
    output:
        netcdf=os.path.join(PROCESSING_DIR, "daily_{year}.nc")
    params:
        year="{year}",
        indir=INDIR,
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        fields=FIELDS
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/%A_%a.out --error=sbatch_dump/%A_%a.err"
    conda:
        "/lustre/soge1/users/spet5107/micromamba/envs/hazGAN-torch"
        # os.path.join("..", "..", CONDA) # or existing named env since 6.14.0 (discouraged)
    log:
        "logs/get_data_{year}.log"
    script:
        "../scripts/get_data.py"
