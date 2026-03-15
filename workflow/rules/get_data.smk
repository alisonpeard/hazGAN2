"""
Get ERA5 data from the SoGE cluster. Run rules from repository root.
"""
def get_param_file():
    """Get a parameter file for the project, or a grid file for the dataset."""
    param_file = RESOURCES_DIR / "params.nc"
    if param_file.exists():
        return param_file
    else:
        raise FileNotFoundError(
            f"Parameter file not found in {str(RESOURCES_DIR)}. "
            "Please create the file first. "
            "This can be empty if no parameters are needed."
        )


rule get_year:
    """
    >>> snakemake --profile profiles/slurm/ --executor slurm --jobs 1 projects/bayofbengal_era5/results/processing/input/2020.nc
    >>> snakemake --profile profiles/cluster --jobs 1 projects/bayofbengal_era5/results/processing/input/2020.nc
    >>> snakemake --profile profiles/cluster --jobs 1 /soge-home/users/spet5107/code/hazGAN2/projects/poweruk_winter/results/processing/input/2005.nc -n
    """
    input:
        indir=INDIR,
        params=get_param_file()
    output:
        netcdf=PROCESSING_DIR / "input" / "{year}.nc"
    params:
        year="{year}",
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        fields=FIELDS,
        dataset=DATASET,
        antecedent_buffer_days=config.get("antecedent_buffer_days", None),
        tmpdir=TMPDIR
    resources:
        cpus_per_task=4,
    conda:
        GEOENV
    log:
        file=Path("logs") / "get" / "{year}.log"
    script:
        Path("..") / "scripts" / "get_data.py"


rule get_all_years:
    """
    Process all input years of input data for the project.

    >>> snakemake --profile profiles/slurm --executor slurm get_all_years --use-conda --jobs 40
    """
    input:
        expand(
            PROCESSING_DIR / "input" / "{year}.nc",
            year=YEARS
        )