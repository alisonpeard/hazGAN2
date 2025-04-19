# hazGAN2 Readme
### Improvements to make / to do
- [ ] Separate conda environments for each rule set
- [ ] Test 2020 subset locally
- [ ] Training on styleGAN2-DA (make very modular)
- [ ] Post-processing samples `process_samples.py`
- [ ] Visualisation
- [ ] Extra analysis code

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
