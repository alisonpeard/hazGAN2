# HazGAN2 readme
This repository contains a `snakemake` workflow to generate multivariate climate data event sets using extreme value theory and generative adversarial networks. The workflow has been made as modular as possible, to facilitate modifications for new applications.

The theory of the workflow is described in [this paper](link.to.paper.com).
### Improvements to make / to do
- [ ] Separate conda environments for each rule set
- [ ] Training on styleGAN2-DA (make very modular)
- [ ] Post-processing samples `process_samples.py`
- [ ] Visualisation rules
- [ ] Extra analysis code/rules
- [ ] Run `get_data` rules locally for 2020 subset
    - [x] `bayofbengal`
    - [ ] `renewablesuk`
- [ ] Run `process_data` rules locally for 2020 subset
    - [ ] `bayofbengal`
    - [ ] `renewablesuk`

### Getting started

Though not always recommended, it is better to submit jobs and set up on a CPU node rather than the head node because head nodes are slow for configuring conda environments etc.
```bash
srun -p Short --pty /bin/bash
```

To set up snakemake on your machine, use the following steps, replacing `conda` with `micromamba` or `mamba` if you prefer:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda config --set channel_priority strict # snakemake complains otherwise
conda install -c conda-forge conda=24.7.1
python -m pip install snakemake-executor-plugin-slurm # snakemake >= 9.0.0, if using SLURM
```

To run a rule, navigate to the repository root, activate snakemake and run the rule. E.g., to 
run all the get_data rules:
```bash
cd hazGAN2
micromamba activate snakemake
snakemake --profile profiles/local/ get_data --use-conda --cores 2
```
or if using SLURM:
```bash
snakemake --profile profiles/cluster/ --executor slurm get_data --use-conda
```

## Modifications and extensions

Touchpoints for modifications:
- `config/config.yaml`: change `project` value
- `config/projects/`: add YAML `{myproject}.yaml` with same structure as existing projects
- `resources/params/`: use `resources/grids/era5.nc` to make `{myproject}.nc` with any spatial parameters for variable construction from raw ERA5 variables (see other param files for examples)
- `workflow/py_utils/era5_utils.py`: add new variable constructions functions defined in `config/projects/{myproject}.yaml`
- `workflow/r_utils/`: add distribution definition files if needed, e.g., `gaussian.R` with functions `cdf()` and `threshold_selector()` following syntax of `genpareto.R` and `weibull.R`

### Data input structure
```
era5dir/
└── {variable_long_name}/
    └── nc/
       └── {variable_long_name}_{year}.nc
 ```

### Repository map
```
hazGAN2/
├── .gitignore
├── README.md
├── docs/
├── profiles/
│   └── {device}/
│       └── config.yaml
├── config/
│   ├── config.yaml
│   └── projects/
│       └── {projectname}.yaml
├── workflow/
│   ├── Snakefile
│   ├── environments/
│   │   ├── {renvs}.yaml
│   │   └── {pythonenvs}.yaml
│   ├── rules/
│   │   └── {rulename}.smk
│   ├── scripts/
│   │   ├── {scriptname}.py
│   │   └── {scriptname}.R
│   ├── py_utils/
│   │   └── {module}.py
│   └── r_utils/
│       └── {module}.R
├── packages/
│   ├── hazGAN/
│   │   ├── pyproject.toml
│   │   └── src/
│   │       └── *
│   └── styleGAN2-DA/
│       ├── pyproject.toml
│       └── src/
│           └── *
├── results/
│   ├── .gitignore
│   └── {projectname}/
│       ├── processing/
│       ├── training/
│       ├── generated/
│       └── analysis/
├── resources/
│   ├── .gitignore
│   ├── grids/
│   │   └── era5.nc
│   └── params/
│       └── {projectname}.nc
│── sbatch_dump/
│   └── .gitignore
└── logs/
    └── .gitignore
```
