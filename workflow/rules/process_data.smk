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
        os.path.join(PROCESSING_DIR, "data_all.nc")

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
        os.path.join("..", "..", CONDA)
    log:
        os.path.join("logs", "resample_{year}.log")
    script:
        os.path.join("..", "scripts", "resample_data.py")


rule concatenate_data:
    input:
        netcdfs=expand(
            os.path.join(PROCESSING_DIR, "resampled_{year}.nc"),
            year=YEARS
        )
    output:
        netcdf=os.path.join(PROCESSING_DIR, "data_all.nc")
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/concat_%A_%a.out --error=sbatch_dump/concat_%A_%a.err"
    conda:
        os.path.join("..", "..", CONDA)
    log:
        os.path.join("logs", "concatenate.log")
    script:
        os.path.join("..", "scripts", "concatenate_data.py")


rule remove_windbombs:
    """Not implementing this."""
    input:
        netcdf=os.path.join(PROCESSING_DIR, "data_all.nc")
    output:
        netcdf=os.path.join(PROCESSING_DIR, "data_nobombs.nc"),
        windbomb=os.path.join(PROCESSING_DIR, "windbomb.npy")
    params:
        threshold=0.82
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/windbombs_%A_%a.out --error=sbatch_dump/windbombs_%A_%a.err"
    conda:
        os.path.join("..", "..", CONDA)
    log:
        os.path.join("logs", "windbombs.log")
    script:
        os.path.join("..", "scripts", "remove_windbombs.py")