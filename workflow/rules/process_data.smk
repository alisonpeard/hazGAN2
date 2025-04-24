"""
Rules to transform daily gridded ERA5 data for RoI into multivariate
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
rule process_all_data:
    """Complete full data processing sequence."""
    input:
        os.path.join(TRAINING_DIR, "jpegs.zip")


checkpoint make_jpegs:
    """Make jpegs of the training data."""
    input:
        data=os.path.join(TRAINING_DIR, "data.nc")
    output:
        outdir=directory(os.path.join(TRAINING_DIR, "jpegs")),
        zipfile=os.path.join(TRAINING_DIR, "jpegs.zip"),
        image_stats=os.path.join(TRAINING_DIR, "image_stats.npz")
    params:
        event_subset=config['event_subset'],
        do_subset=False, # TODO: do this better
        eps = 1e-6,
        domain = config["domain"],
        resx = RESOLUTION['lon'],
        resy = RESOLUTION['lat'],
    conda:
        PYENV
    log:
        file=os.path.join("logs", "make_jpegs.log")
    script:
        os.path.join("..", "scripts", "make_jpegs.py")


rule make_training_data:
    """Convert dataframe to training netCDF."""
    input:
        data_all=os.path.join(PROCESSING_DIR, "data_all.nc"),
        events=os.path.join(PROCESSING_DIR, "events.parquet"),
        metadata=os.path.join(PROCESSING_DIR, "event_metadata.csv"),
        # medians=os.path.join(PROCESSING_DIR, "medians.csv")
    output:
        data=os.path.join(TRAINING_DIR, "data.nc")
    params:
        fields=FIELDS
    conda:
        PYENV
    log:
        file=os.path.join("logs", "make_training.log")
    script:
        os.path.join("..", "scripts", "make_training.py")


rule fit_marginals:
    """Fit semi-parametric marginals to the data along the time dimension.
    
    Usage:
    >>> snakemake --profile profiles/local/ fit_marginals --use-conda --cores 2
    """
    input:
        metadata=os.path.join(PROCESSING_DIR, "event_metadata.csv"),
        daily=os.path.join(PROCESSING_DIR, "daily.parquet")
    output:
        events=os.path.join(PROCESSING_DIR, "events.parquet")
    params:
        fields=FIELDS,
        q=MTHRESH
    conda:
        RENV
    log:
        file="logs/fit_marginals.log"
    script:
        os.path.join("..", "scripts", "fit_marginals.R")


rule extract_events:
    """Remove seasonality and extract storm events from the data.

    TODO: may need to re-add medians later. Calculating medians on
    the fly is less robust.
    
    Params:
        sfunc: str
            Any function that takes args df and vars (list).
            Must be defined in utils.R.
        rfunc: str
            Any function that can be fed into aggregate(). Standard
            or defined in utils.R.
    """
    input:
        # os.path.join(".snakemake", "conda", ".rpot_installed"),
        netcdf=os.path.join(PROCESSING_DIR, "data_all.nc")
    output:
        metadata=os.path.join(PROCESSING_DIR, "event_metadata.csv"),
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
        file=os.path.join("logs", "extract_events.log")
    script:
        os.path.join("..", "scripts", "extract_events.R")


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
        file=os.path.join("logs", "concatenate.log")
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
        file=os.path.join("logs", "resample_{year}.log")
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