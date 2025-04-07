## Environment set up
To set up snakemake on your machine you can use the following steps
```bash
micromamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge conda=24.7.1
install snakemake-executor-plugin-slurm
```

Then to run a rule, navigate to the root of the repository and run:
```bash
micromamba activate snakemake
# your rule
```

## Getting the data
Data needs to collected from the SoGE cluster using the rules in `get_data.smk` rule. ERA5 data is stored in `/soge-home/data/analysis/era5/0.28125x0.28125/hourly/`. If you don't have access to the SoGE cluster or you don't want to use ERA5 you will need to modify this rule and its corresponding script `get_data.py`.

Project-level settings for processing the input data are defined in `config/projects/<project>.yaml`. The fields section defines the three climate fields that will be used in the event footprints. These are defined along with their `args` (raw ERA5 fields to construct them from), `func` (function to construct them), and `agg` function for their temporal aggregation ($h_{k|t}(x), k=1,2,3$ in the paper). The function must be defined by name in `workflow/scripts/era5_utils.py`. Currently, `agg` only accepts standard aggregation variables but this can be changed.

To run this rule (on SLURM), follow the setup steps above and run:
```bash
snakemake --profile profiles/cluster/ --executor slurm get_data --use-conda
```