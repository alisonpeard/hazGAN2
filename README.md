# HazGAN2 readme
This repository contains a `snakemake` workflow to generate multivariate climate event sets using extreme value theory and generative adversarial networks.

The workflow has been made as modular as possible to facilitate modifications for new applications.

The theory of the workflow is described in [this paper](link.to.paper.com) and the rest of this readme describes how to get started with the workflow.

## Getting started

```bash
# clone the repository
git clone git@github.com:alisonpeard/hazGAN2.git
cd hazGAN2
```

It is better to run snakemake from a CPU node rather than the head node, head nodes are extremely slow for creating conda environments.
```bash
# login to a CPU node
srun -p Short --pty /bin/bash
```

To set up snakemake on your machine, enter the following in the terminal, replacing `conda` with `micromamba` or `mamba` if you prefer:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda config --set channel_priority strict # snakemake complains otherwise
conda install -c conda-forge conda=24.7.1
python -m pip install snakemake-executor-plugin-slurm # snakemake >= 9.0.0, if using SLURM
```

### Running rules
To run a rule, navigate to the repository root, activate snakemake and run the rule. E.g., to 
run all the rules in `get_data.smk` locally:
```bash
cd hazGAN2
micromamba activate snakemake
snakemake --profile profiles/local/ get_all_data --use-conda --cores 2
```
or if using SLURM:
```bash
snakemake --profile profiles/cluster/ --executor slurm get_all_data --use-conda
```
#### Sample rules:
```bash
snakemake --profile profiles/local/ process_all_data --use-conda --cores 2
snakemake --profile profiles/local/ fit_marginals --use-conda --cores 2
```
> 💭 For Apple Silicon, the R package `r-extremes` is not available on the conda `osx-arm64` subdirectory, so installation must be manually set to the `osx-64` subdirectory. If running a rule that will install the R environment, prefix the command with `CONDA_SUBDIR=osx-64`, e.g.,
```bash
CONDA_SUBDIR=osx-64 snakemake --profile profiles/local/ fit_marginals --use-conda --cores 2
```

### Making reports and DAGs
and to output the DAG for a specific rule:
```bash
# without files
snakemake process_all_data --dag | dot -Tpdf > docs/process_all_data.pdf
snakemake process_all_data --dag | dot -Tsvg > docs/process_all_data.svg

# with files
snakemake process_all_data --filegraph | dot -Tpdf > docs/process_all_data.pdf
snakemake process_all_data --filegraph | dot -Tsvg > docs/process_all_data.svg

# report
snakemake process_all_data --report docs/process_all_data.html
```
## Modifications and extensions
#### Basic
To create a new project you need to make the following changes:
1. `config/config.yaml`: change `project` value
2. `config/projects/`: add a YAML file named `{myproject}.yaml` with the same structure as existing project YAMLs
3. `resources/params/`: use `resources/grids/era5.nc` to make `{myproject}.nc` in `resources/params/` with any spatial parameters for variable construction from raw ERA5 variables (see other param files for examples).

#### Advanced
The following files may also need to be modified, depending on the project requirements:
- Add new variable construction functions to
    - `workflow/scripts/era5_utils.py`
    - `config/projects/{myproject}.yaml`
- Add new temporal aggregation functions $h_{k|t}(x)$ to
    - `workflow/scripts/era5_utils.py`
    - `config/projects/{myproject}.yaml:hfunc`
- Add new deasonalization $s_{|t}(x)$ methods to
    - `worflow/r_utils/sfuncs.R`
    - `config/projects/{myproject}.yaml:sfunc`
- Add new risk functionals $r_{|ijk}(x)$ to
    - `workflow/r_utils/r)funcs.R`
    - `config/projects/{myproject}.yaml:rfunc`
- Add new distribution definitions to 
    - `workflow/r_utils/{distribution}.R`:
        - `cdf()`
        - `threshold_selector()`


## Further information
### Data input structure
If you don't have access to the SoGE filestore, you should set up the input data structure as follows:
```
input/
└── {variable_long_name}/
    └── nc/
       └── {variable_long_name}_{year}.nc
 ```
 and the output (`results`) folder has the following structure:
 ```bash
 :
 ├── results/
 │   ├── .gitignore
 │   └── {projectname}/
 │       ├── processing/
 │       ├── training/
 │       ├── generated/
 │       └── analysis/
 :
```

### Repository map
The full repository is structured as follows:
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

### Tasks
- [ ] Separate conda environments for each rule set
- [ ] Training on styleGAN2-DA (make very modular)
- [ ] Post-processing samples `process_samples.py`
- [ ] Visualisation rules
- [ ] Extra analysis code/rules
- [ ] Either add test set back in or full remove
- [ ] Unit tests
- [ ] Create fallback `params` file
- [ ] Run `get_data` rules locally for 2020 subset
    - [x] `bayofbengal`
    - [ ] `renewablesuk`
- [ ] Run `process_data` rules locally for 2020 subset
    - [ ] `bayofbengal`
    - [ ] `renewablesuk`