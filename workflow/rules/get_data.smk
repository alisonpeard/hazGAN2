"""
Get ERA5 data from SoGE cluster. Run from repository root.

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
rule get_all_data:
    """Rule to process all years for the project."""
    input:
        expand(
            os.path.join(PROCESSING_DIR, "daily_{year}.nc"),
            year=YEARS
        )

rule get_year:
    input:
        indir=INDIR,
        params=os.path.join(RESOURCES_DIR, "params", f"{PROJECT}.nc")
    output:
        netcdf=os.path.join(PROCESSING_DIR, "daily_{year}.nc")
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
        PYENV
    log:
        file="logs/get_{year}.log"
    script:
        "../scripts/get_data.py"
