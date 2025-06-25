"""
Rules to transform daily gridded ERA5 data for RoI into multivariate
event footprints for training GAN.

# environment set up
>>> micromamba create -c conda-forge -c bioconda -n snakemake snakemake
>>> conda activate snakemake
>>> conda install -c conda-forge conda=24.7.1
>>> python -m pip install snakemake-executor-plugin-slurm
 
# run rule project root and login node
>>> micromamba activate snakemake
>>> snakemake --profile profiles/cluster/ --executor slurm process_data --use-conda
"""
rule process_all_data:
    """Complete full data processing sequence."""
    input:
        os.path.join(TRAINING_DIR, "images.zip")


checkpoint make_rgb_images:
    """Make PNGs of the training data."""
    input:
        data=os.path.join(TRAINING_DIR, "data.nc")
    output:
        outdir=directory(os.path.join(TRAINING_DIR, "rgb")),
        zipfile=os.path.join(TRAINING_DIR, "images.zip"),
        image_stats=os.path.join(TRAINING_DIR, "image_stats.npz")
    params:
        event_subset=config['event_subset'],
        do_subset=True, # TODO: do this better
        eps = 1e-6,
        domain = config["domain"],
        resx = RESOLUTION['lon'],
        resy = RESOLUTION['lat'],
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "make_rgb.log")
    script:
        os.path.join("..", "scripts", "make_rgb_images.py")


rule make_training_data:
    """Convert dataframe to training netCDF."""
    input:
        events=os.path.join(PROCESSING_DIR, "events.parquet"),
        metadata=os.path.join(PROCESSING_DIR, "event_metadata.parquet"),
        medians=os.path.join(PROCESSING_DIR, "medians.parquet")
    output:
        data=os.path.join(TRAINING_DIR, "data.nc")
    params:
        fields=FIELDS
    conda:
        GEOENV
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
        metadata=os.path.join(PROCESSING_DIR, "event_metadata.parquet"),
        daily=os.path.join(PROCESSING_DIR, "daily.parquet")
    output:
        events=os.path.join(PROCESSING_DIR, "events.parquet")
    params:
        fields=FIELDS
    conda:
        RENV
    log:
        file="logs/fit_marginals.log",
        level="DEBUG"
    script:
        os.path.join("..", "scripts", "fit_marginals.R")


rule extract_events:
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
        medians=os.path.join(PROCESSING_DIR, "medians.parquet"),
        metadata=os.path.join(PROCESSING_DIR, "event_metadata.parquet"),
        daily=os.path.join(PROCESSING_DIR, "daily.parquet")
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
        cpus_per_task=4
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
            os.path.join(PROCESSING_DIR, "resampled", "{year}.nc"),
            year=YEARS
        )
    output:
        netcdf=os.path.join(PROCESSING_DIR, "data_all.nc")
    params:
        exclude=config['exclude']
    resources:
        cpus_per_task=4
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "concatenate.log")
    script:
        os.path.join("..", "scripts", "concatenate_data.py")


rule resample_year:
    """Resample the data to the desired resolution."""
    input:
        netcdf=os.path.join(PROCESSING_DIR, "input", "{year}.nc")
    output:
        netcdf=os.path.join(PROCESSING_DIR, "resampled", "{year}.nc")
    params:
        year="{year}",
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
        fields=FIELDS
    conda:
        GEOENV
    resources:
        cpus_per_task=4
    log:
        file=os.path.join("logs", "resample", "{year}.log")
    script:
        os.path.join("..", "scripts", "resample_data.py")