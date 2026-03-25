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
>>> snakemake --profile profiles/cluster/ --use-conda process_all_data
>>> snakemake --profile profiles/slurm/ --executor slurm --use-conda process_all_data
"""
from pathlib import Path


rule resample_year:
    """Resample the data to the desired resolution.
    >>> snakemake --profile profiles/cluster/ projects/poweruk/results/processing/resampled/2017.nc
    """
    input:
        netcdf=PROCESSING_DIR / "input" / "{year}.nc"
    output:
        netcdf=PROCESSING_DIR / "resampled" / "{year}.nc"
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
        file=Path("logs") / "resample" / "{year}.log"
    script:
        Path("..") / "scripts" / "resample_data.py"


rule concatenate_resampled:
    """Concatenate all the years into a single netcdf file."""
    input:
        netcdfs=expand(
            PROCESSING_DIR / "resampled" / "{year}.nc",
            year=YEARS
        )
    output:
        netcdf=PROCESSING_DIR / "resampled_all.nc"
    params:
        exclude=EXCLUDE
    resources:
        cpus_per_task=4
    conda:
        GEOENV
    log:
        file=Path("logs") / "concatenate.log"
    script:
        Path("..") / "scripts" / "concatenate_data.py"



rule extract_events:
    """Remove seasonality and extract storm events from the data.
    
    Params:
        sfunc: str
            Any function that takes args df and vars (list).
            Must be defined in utils.R.
        rfunc: str
            Any function that can be fed into aggregate(). Standard
            or defined in utils.R.
        R_funcs: 
    
    >>> snakemake --profile profiles/cluster projects/poweruk2/results/processing/events.parquet
    """
    input:
        netcdf=PROCESSING_DIR / "resampled_all.nc"
    output:
        medians=PROCESSING_DIR / "climatology.parquet",
        metadata=PROCESSING_DIR / "event_metadata.parquet",
        daily=PROCESSING_DIR / "event_cubes.parquet"
    params:
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
        xmin=LONGITUDE['min'],
        xmax=LONGITUDE['max'],
        ymin=LATITUDE['min'],
        ymax=LATITUDE['max'],
        fields=FIELDS,
        rfunc=RFUNC,
        sfunc=SFUNC,
        R_funcs=RFUNCS
    resources:
        cpus_per_task=4
    conda:
        RENV
    log:
        file=Path("logs") / "extract_events.log"
    script:
        Path("..") / "scripts" / "extract_events.R"


rule fit_marginals:
    """Fit semi-parametric marginals to the data along the time dimension.
    
    Usage:
    >>> snakemake --profile profiles/local/ fit_marginals --use-conda --cores 2
    >>> snakemake --profile profiles/cluster projects/poweruk2/results/processing/events.parquet
    """
    input:
        metadata=PROCESSING_DIR / "event_metadata.parquet",
        cubes=PROCESSING_DIR / "event_cubes.parquet"
    output:
        events=PROCESSING_DIR / "event_footprints.parquet"
    params:
        fields=FIELDS,
        R_funcs=RFUNCS
    conda:
        RENV
    log:
        file=Path("logs") / "fit_marginals.log",
        level="DEBUG"
    script:
        Path("..") / "scripts" / "fit_marginals.R"



rule make_training_data:
    """Convert dataframe to training netCDF.
    
    >>> snakemake --profile profiles/slurm make_training_data
    """
    input:
        events=PROCESSING_DIR / "event_footprints.parquet",
        metadata=PROCESSING_DIR / "event_metadata.parquet",
        medians=PROCESSING_DIR / "climatology.parquet"
    output:
        data=PROCESSING_DIR / "data.nc"
    params:
        fields=FIELDS,
        domain=config["domain"]
    conda:
        GEOENV
    log:
        file=Path("logs") / "make_training.log"
    script:
        Path("..") / "scripts" / "make_training.py"


checkpoint make_rgb_images:
    """Make PNGs of the training data.
    
    >>> snakemake --profile profiles/slurm make_rgb_images
    """
    input:
        data=PROCESSING_DIR / "data.nc"
    output:
        outdir=directory(Path(TRAINING_DIR) / "rgb"),
        zipfile=Path(TRAINING_DIR) / "images.zip",
        image_stats=Path(TRAINING_DIR) / "image_stats.npz"
    params:
        subset=config['event_subset'],
        rpmax=1e6,
        eps=1e-6,
        domain=config["domain"],
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
    conda:
        GEOENV
    log:
        file=Path("logs") / "make_rgb.log"
    script:
        Path("..") / "scripts" / "make_rgb_images.py"


rule process_all_data:
    """Complete full data processing sequence."""
    input:
        Path(TRAINING_DIR) / "images.zip"
