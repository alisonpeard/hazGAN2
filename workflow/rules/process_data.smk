"""
Get data from SoGE cluster. Run from repository root.

```bash
# environment set up
micromamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge conda=24.7.1
python -m pip install snakemake-executor-plugin-slurm

# run rule project root and login node
micromamba activate snakemake
snakemake --profile profiles/cluster/ --executor slurm process_data --use-conda
```
"""
rule process_data:
    """Rule to process all years for the project."""
    input:
        expand(
            os.path.join(PROCESSING_DIR, "resampled_{year}.nc"),
            year=YEARS
        )

rule resample_year:
    input:
        netcdf=os.path.join(PROCESSING_DIR, "daily_{year}.nc")  
    output:
        netcdf=os.path.join(PROCESSING_DIR, "resampled_{year}.nc")
    params:
        year="{year}",
        resx=config["resolution"]['lon'],
        resy=config["resolution"]['lat'],
        fields=FIELDS
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/resample_%A_%a.out --error=sbatch_dump/resample_%A_%a.err"
    conda:
        # "/lustre/soge1/users/spet5107/micromamba/envs/hazGAN-torch"
        os.path.join("..", "..", CONDA)
    log:
        "logs/resample_{year}.log"
    script:
        "../scripts/resample_data.py"