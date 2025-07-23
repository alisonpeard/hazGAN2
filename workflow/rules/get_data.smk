"""
Get ERA5 data from SoGE cluster. Run from repository root.

NOTE: add GCS bucket methods too.

# environment set up
>>> micromamba create -c conda-forge -c bioconda -n snakemake snakemake
>>> conda activate snakemake
>>> >>> conda install -c conda-forge conda=24.7.1
install snakemake-executor-plugin-slurm

# run rule from project root while on the login node
>>> micromamba activate snakemake
>>> snakemake --profile profiles/slurm/ --executor slurm get_data --use-conda --jobs 20
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
    """Rule to process all years for the project.
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
    """
    input:
        indir=INDIR,
        params=get_param_file()
    output:
        netcdf=os.path.join(PROCESSING_DIR, "input", "{year}.nc")
    params:
        year="{year}",
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        timecol=TIMECOL,
        fields=FIELDS,
        dataset=DATASET
    resources:
        cpus_per_task=4,
        # slurm_extra="--output=sbatch_dump/get_%A_%a.out --error=sbatch_dump/get_%A_%a.err" 
    conda:
        GEOENV
    log:
        file="logs/get/{year}.log"
    script:
        "../scripts/get_data.py"
