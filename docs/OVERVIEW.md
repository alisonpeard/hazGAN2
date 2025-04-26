## Environment set up

Though not always recommended, it is better to submit jobs and set up on a CPU node rather than the head node because head nodes are slow for configuring conda environments etc.:
```bash
srun -p Short --pty /bin/bash
```

To set up snakemake on your machine, use the following steps:
```bash
micromamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda config --set channel_priority strict
conda install -c conda-forge conda=24.7.1
python -m pip install snakemake-executor-plugin-slurm # snakemake >= 9.0.0
```

Then to run a rule, navigate to the root of the repository (`hazGAN`) and run:
```bash
micromamba activate snakemake
# your rule
```

To view an overview of the DAG for a rule:
```bash
snakemake --profile profiles/cluster/ process_data --dag dot | dot -Tpng > workflow/dag.png
```

> ‼️ There's currently an issue using existing conda environments for more than one rule. It's easier to set the `conda` key to a path to a YAML file.

## Getting the data
### Data source
Data is collected from the SoGE cluster using rules defined in `get_data.smk`.
* ERA45 data is stored in `/soge-home/data/analysis/era5/0.28125x0.28125/hourly/`
* Other data is stored in `/soge-home/data/` but rules have not been defined for this yet.
* If you don't have access to the SoGE cluster or you don't want to use ERA5, you will need to modify the `get_data.smk` rule and its companion script `get_data.py`.

### Data formatting
Each project has its own data formatting rules. These are defined in `config/projects/<project name>.yaml`. The following fields are required:
* `year0`, `yearn` : start and end years for the project
* `longitude` : longitude bounds for the project with keys `min` and `max`
* `latitude` : latitude bounds for the project with keys `min` and `max`
* `fields`: three climate fields for the project. Each key has the following sub-keys:
  * `args`: list of raw ERA5 fields to construct the field from, defined in `workflow/scripts/era5_utils.py`
  * `func`: function to apply to the arguments, must be defined in `workflow/scripts/era5_utils.py`
  * `agg`: aggregation function to apply to the output of the function, must be standard aggregation variables (e.g. `mean`, `max`, `min`, etc.)

To run this rule (on SLURM):
1. Follow setup instructions above
2. Navigate to the root of the repository (`hazGAN`)
3. Run the following command:
    ```bash
    snakemake --profile profiles/cluster/ --executor slurm get_data --use-conda
    ```

> 💡 The first time this is run, snakemake will create a conda environment for each rule. This can take a while. Once the environment is created, it will be cached and reused for subsequent runs.


## Processing the data
### 1. Resampling to desired resolution
```bash
snakemake --profile profiles/cluster/ --executor slurm get_data --use-conda
```
> ‼️ This didn't re-run `get_data` when I changed the fields in `config/projects/<project name>.yaml`.

### 2. Processing resampled data for analysis

This consists of two rules:
- `concatenate_data`: concatenates the data from all years into a single file and assignes grid indices
- `remove_windbombs`: specific to `bayofbengal` project, this removes windbombs from the data.

### 3. Extracting event footprints

This rules requires Rscripts.

### 4. Fitting marginal distributions

This rules requires Rscripts.

### 5. Making netCDFs for custom training
### 6. Making JPEGs for StyleGAN2
