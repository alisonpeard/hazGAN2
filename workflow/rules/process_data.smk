"""
Rules to turn daily gridded ERA5 data for RoI into multivariate
event footprints for training GAN.

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
    """Complete full data processing sequence."""
    input:
        os.path.join(PROCESSING_DIR, "storms.parquet")


checkpoint make_jpegs:
    """Make jpegs of the training data."""
    input:
        data=os.path.join(TRAINING_DIR, "data.nc")
    output:
        outdir=os.path.join(TRAINING_DIR, "jpegs"),
        image_stats=os.path.join(TRAINING_DIR, "image_stats.npz")
    params:
        threshold=config['event_threshold'],
        eps = 1e-6,
        domain = "gumbel",
        resx = RESOLUTION['lon'],
        resy = RESOLUTION['lat'],
    conda:
        PYENV
    script:
        os.path.join("..", "scripts", "make_jpegs.py")

rule make_training_data:
    """Make training data for GAN."""
    input:
        data_all=os.path.join(PROCESSING_DIR, "data_all.nc"),
        storms=os.path.join(PROCESSING_DIR, "storms.parquet"),
        metadata=os.path.join(PROCESSING_DIR, "storm_metadata.csv"),
        medians=os.path.join(PROCESSING_DIR, "medians.csv")
    output:
        data=os.path.join(TRAINING_DIR, "data.nc")
    params:
        fields=FIELDS
    conda:
        PYENV
    script:
        os.path.join("..", "scripts", "make_training.py")

rule fit_marginals:
    """Fit semi-parametric marginals to the data along the time dimension."""
    input:
        medians=os.path.join(PROCESSING_DIR, "medians.csv"),
        metadata=os.path.join(PROCESSING_DIR, "storm_metadata.csv"),
        daily=os.path.join(PROCESSING_DIR, "daily.parquet"),
    output:
        storms=os.path.join(PROCESSING_DIR, "storms.parquet"),
    params:
        fields=FIELDS,
        q=0.95
    conda:
        RENV
    script:
        os.path.join("..", "scripts", "fit_marginals.R")


rule extract_storms:
    """Remove seasonality and extract storm events from the data.

    Params:
        sfunc: str
            Any function that takes args df and vars (list).
            Must be defined in utils.R.
        rfunc: str
            Any function that can be fed into aggregate(). Standard
            or defined in utils.R.
    """
    input:
        netcdf=os.path.join(PROCESSING_DIR, "data_all.nc")
    output:
        medians=os.path.join(PROCESSING_DIR, "medians.csv"),
        metadata=os.path.join(PROCESSING_DIR, "storm_metadata.csv"),
        daily=os.path.join(PROCESSING_DIR, "daily.parquet"),
    params:
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
        xmin=LONGITUDE['min'],
        xmax=LONGITUDE['max'],
        ymin=LATITUDE['min'],
        ymax=LATITUDE['max'],
        fields=FIELDS,
        rfunc=RFUNC,
        sfunc=SFUNC
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/storms_%A_%a.out --error=sbatch_dump/storms_%A_%a.err"
    conda:
        RENV
    log:
        os.path.join("logs", "extract_storms.log")
    script:
        os.path.join("..", "scripts", "extract_storms.R")

rule concatenate_data:
    """Concatenate all the years into a single netcdf file."""
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
        PYENV
    log:
        os.path.join("logs", "concatenate.log")
    script:
        os.path.join("..", "scripts", "concatenate_data.py")


rule resample_year:
    """Resample the data to the desired resolution."""
    input:
        netcdf=os.path.join(PROCESSING_DIR, "daily_{year}.nc")  
    output:
        netcdf=os.path.join(PROCESSING_DIR, "resampled_{year}.nc")
    params:
        year="{year}",
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
        fields=FIELDS
    conda:
        PYENV
    resources:
        cpus_per_task=4,
        slurm_extra="--output=sbatch_dump/resample_%A_%a.out --error=sbatch_dump/resample_%A_%a.err"
    log:
        os.path.join("logs", "resample_{year}.log")
    script:
        os.path.join("..", "scripts", "resample_data.py")

# rule remove_windbombs:
#     """Not implemented."""
#     input:
#         netcdf=os.path.join(PROCESSING_DIR, "data_all.nc")
#     output:
#         netcdf=os.path.join(PROCESSING_DIR, "data_nobombs.nc"),
#         windbomb=os.path.join(PROCESSING_DIR, "windbomb.npy")
#     params:
#         threshold=0.82
#     resources:
#         cpus_per_task=4,
#         slurm_extra="--output=sbatch_dump/windbombs_%A_%a.out --error=sbatch_dump/windbombs_%A_%a.err"
#     conda:
#         os.path.join("..", "..", CONDA)
#     log:
#         os.path.join("logs", "windbombs.log")
#     script:
#         os.path.join("..", "scripts", "remove_windbombs.py")