"""
Get ERA5 data from the SoGE cluster. Run rules from repository root.
"""
def get_param_file():
    """Get a parameter file for the project, or a grid file for the dataset."""
    if os.path.exists(os.path.join(RESOURCES_DIR, "params.nc")):
        return os.path.join(RESOURCES_DIR, "params.nc")
    else:
        raise FileNotFoundError(
            f"Parameter file not found in {RESOURCES_DIR}. "
            "Please create the file first. This can be empty if no paramaters are needed."
        )


rule get_all_years:
    """
    Process all input years of input data for the project.

    >>> snakemake --profile profiles/slurm/ --executor slurm get_all_years --use-conda --jobs 40
    """
    input:
        expand(
            os.path.join(PROCESSING_DIR, "input", "{year}.nc"),
            year=YEARS
        )


rule get_year:
    """
    >>> snakemake --profile profiles/slurm/ --executor slurm --jobs 1 projects/bayofbengal_era5/results/processing/input/2020.nc
    >>> snakemake --profile profiles/cluster --jobs 1 projects/bayofbengal_era5/results/processing/input/2020.nc
    >>> snakemake --profile profiles/cluster --jobs 1 projects/poweruk/results/processing/input/2017.nc
    """
    input:
        indir=INDIR,
        params=get_param_file()
    output:
        netcdf=os.path.join(PROCESSIccleNG_DIR, "input", "{year}.nc")
    params:
        year="{year}",
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        timecol=TIMECOL,
        fields=FIELDS,
        dataset=DATASET,
        antecedent_buffer_days=config.get("antecedent_buffer_days", None)
    resources:
        cpus_per_task=4,
        # slurm_extra="--output=sbatch_dump/get_%A_%a.out --error=sbatch_dump/get_%A_%a.err" 
    conda:
        GEOENV
    log:
        file="logs/get/{year}.log"
    script:
        "../scripts/get_data.py"