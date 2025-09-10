[![Snakemake](https://img.shields.io/badge/snakemake->=8.0.0-9cf.svg?style=flat)](https://snakemake.readthedocs.io)

[![Python](https://img.shields.io/badge/python-3.12.9-9cf.svg?style=flat)](https://snakemake.readthedocs.io)

# HazGAN2 readme

This repository contains a [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to generate multivariate climate event sets using extreme value theory and generative adversarial networks.

The workflow has been made as modular as possible to facilitate modifications for new applications.

StyleGAN notes: https://github.com/NVlabs/stylegan2-ada-pytorch/issues/11

The theory of the workflow is described in [this paper](link/to/paper.com) and the rest of this README describes basic use.

## Current status [keep updated]

Date: 08-07-2025

- **Complete:**

    - Finished Bay of Bengal workflow ✅

- **In progress:**

    - Modifying workflow to take in event dates

> ### 💡 Note on snakemake relative paths

> * Input, output, log, and benchmark files are considered to be relative to the working directory (either the directory in which you have invoked Snakemake or whatever was specified for --directory or the workdir: directive).
> * Any other directives (e.g. conda:, include:, script:, notebook:) consider paths to be relative to the Snakefile they are defined in.

## Quick start

```bash
micromamba create -c conda-forge -c bioconda -n snakemake snakemake conda=24.7.1 -y
micromamba activate snakemake
# conda config --set-channel_priority strict  # (or change in .condarc file)
python -m pip install snakemake-executor-plugin-slurm # snakemake >= 9.0.0, if using 

# set up the required conda environments
snakemake --profile profiles/cluster --conda-create-envs-only

# Run the entire workflow
snakemake --profile profiles/local -n
snakemake --profile profiles/local
```

## Things to do outside of this workflow

This workflow is just for generating event sets, to keep it clean and modular, any downstream analysis should be done externally with code in the `project/analysis` directory. This includes:

- Generating a `params.nc` file for the project
- All extra analysis (in scripts/Jupyter notebooks)


### Running rules

To run a rule on your local machine, navigate to the repository root, activate snakemake and run the rule. E.g., to 
run all the rules in `get_data.smk` locally:

```bash
cd hazGAN2
micromamba activate snakemake
snakemake --profile profiles/local get_all_data
```

You may need to modify `config/config.yaml` and to point to the correct input data location and Python environment definition files (in `workflow/environments`).

#### Using R environments on Apple Silicon

For Apple Silicon, the R package `r-extremes` is not available on the conda `osx-arm64` subdirectory, so installation must be manually set to the `osx-64` subdirectory. If running a rule that will install the R environment for the first time, prefix the command with `CONDA_SUBDIR=osx-64`, e.g.,
```bash
CONDA_SUBDIR=osx-64 snakemake --profile profiles/local/ process_all_data --use-conda --cores 2
```

#### Running rules on the cluster

When running for the first time, it is better to run snakemake from a CPU node rather than the head node, head nodes are extremely slow for creating conda environments. After that, you can run snakemake from the head node (required for SLURM). To do this, use the following command to login to a CPU node:

```bash
snakemake --profile profiles/cluster/ my_rule --conda-create-envs-only
```

The `config.yaml` in `profiles/local/` contains snakemake command line arguments that are usually run when using the local machine. There are also profiles for running on the cluster in `profiles/cluster/`, and running on the cluster with SLURM in `profiles/slurm`. To run the rule on the cluster use the following command:

```bash
snakemake --profile profiles/cluster/ get_all_data
```
or to send it to SLURM
```bash
snakemake --profile profiles/slurm/ get_all_data
```
You can modify any of the profiles or make a new one to suit your needs.

> 💭 The SoGE cluster has pre-defined SLURM defaults per user and doesn't allow users to set their account. The `snakemake-executor-plugin-slurm` will always attempt to set an account. To override this behavior, the `--slurm-no-account` flag should be used.

Also, you may need a snakemake process to run for a while while you have a job running, you can sbatch the manager script to run in the background:
```bash
micromamba activate snakemake
cd mistral/alison/hazGAN2
sbatch --job-name=snakemake_manager --output=sbatch_dump/snakemake_manager_%j.out --error=sbatch_dump/snakemake_manager_%j.err --wrap="snakemake --profile profiles/slurm/ process_all_data
```
#### Sample rules:
```bash
snakemake --profile profiles/local/ process_all_data --use-conda --cores 2
snakemake --profile profiles/local/ fit_marginals --use-conda --cores 2
```

It's better not to use a head node on the cluster as they are really slow, snakemake will give this warning—-ignore it:
```bash
You are running snakemake in a SLURM job context. This is not recommended, as it may lead to unexpected behavior. Please run Snakemake directly on the login node.
```

### Making reports and DAGs
and to output the DAG for a specific rule:
```bash
# without files
snakemake process_all_data --dag | dot -Tpdf > docs/process_all_data.pdf
snakemake process_all_data --dag | dot -Tsvg > docs/process_all_data.svg

snakemake --profile profiles/cluster all --dag | dot -Tpdf > docs/workflow.pdf

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
    
        Note that the parameters will always be named (`thresh`, `scale`, `shape`), regardless of the distribution. In cases where this doesn't match, assign each parameters to the most appropriate names and set the other to `NA`, e.g., for a normal distribution `thresh` is the mean, `scale` is the standard deviation, and `shape` is not used.


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
│       │   └── yearly/
│       ├── training/
│       │   └── rgb/
│       ├── generated/
│       │   └── stylegan-output/
│       │       └── {networkname}/
│       │               └── generated-{domain}/
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
### DVC data versioning

To set up DVC, install it into its own conda environment:
```bash
# Create new environment
conda create -n dvc python=3.10
conda activate dvc

# Install OpenSSL through conda (instead of relying on system OpenSSL)
conda install -c conda-forge openssl

pip install "dvc[gdrive]"
```

If the above doesn't work, try installing DVC through conda-forge instead:
```bash
conda install -c conda-forge dvc
pip install "dvc[gdrive]"  # To add Google Drive support
```

The original data for [the paper](link/to/paper) is versioned [here](https://github.com/alisonpeard/hazGAN-data) and can be cloned independentendly using

```bash
git clone git@github.com:alisonpeard/hazGAN-data.git
dvc pull
```

You will need to sign into your Google account and grant access to the DVC app. If you are SSH-ing into a remote server, you need to
1. Connect with local port forwarding
```bash
ssh -L 8080:localhost:8080 <username>@<remote_server>
```
2. After running `dvc pull`, copy and paste the URL from the terminal into your local browser. This will open a Google login page. After logging in, you will be redirected to a page with a code. Copy this code and paste it into the terminal on the remote server.
3. After pasting the code, you will be asked to grant access to the DVC app. Click "Allow" and you should see a message saying "You can now close this window". You can now close the browser window and return to the terminal.
4. You should now see a message saying "Authentication successful".

To set up your own DVC file tracking, in the parent repo type `dvc init` and `dvc add results`, to start tracking everything in the results folder with DVC. DVC will prompt you to commit these changes to git. Follow the instructions [here](https://dvc.org/doc/user-guide/data-management/remote-storage/google-drive#using-a-custom-google-cloud-project-recommended) to set up a Google Cloud Project and link it to your DVC repo. You need to set the project status to published rather than testing to allow access. Set the upstream remote and push.

> 💭 There is also an `ssh` version of this for remote filestores. 


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
    - [ ] `poweruk`uk`
