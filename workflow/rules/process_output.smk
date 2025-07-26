rule process_generated:
    """Transform generated images to netCDF and apply inverse
    transformations.
    
    >>> snakemake --profile profiles/slurm process_generated
    """
    input:
        image_dir=os.path.join(GENERATED_DIR, "images"),
        image_stats=os.path.join(TRAINING_DIR, "image_stats.npz"),
        training_dir=os.path.join(TRAINING_DIR, "rgb"),
        training_data=os.path.join(TRAINING_DIR, "data.nc")
    output:
        netcdf=os.path.join(GENERATED_DIR, "netcdf", "data.nc"),
        train=os.path.join(GENERATED_DIR, "netcdf", "train.nc")
    params:
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
        do_subset=True,
        event_subset=config['event_subset'],
        fields=FIELDS,
        domain=config["domain"]
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "process_generated.log")
    script:
        os.path.join("..", "scripts", "process_generated.py")
        

rule make_benchmarks:
    """Create benchmark datasets with assumption of total independence/dependence.

    NOTE: Sampling from base distribution of events (not extremes).
    """
    input:
        data=os.path.join(TRAINING_DIR, "data.nc")
    output:
        dependent=os.path.join(GENERATED_DIR, "netcdf", "dependent.nc"),
        independent=os.path.join(GENERATED_DIR, "netcdf", "independent.nc")
    params:
        resx=RESOLUTION['lon'],
        resy=RESOLUTION['lat'],
        year0=YEAR0,
        yearn=YEARN,
        nyrs=500,
        fields=FIELDS
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "make_benchmarks.log")
    script:
        os.path.join("..", "scripts", "make_benchmarks.py")
