[![Python](https://img.shields.io/badge/python-3.12.9-9cf.svg?style=flat)](https://snakemake.readthedocs.io) [![Snakemake](https://img.shields.io/badge/snakemake->=8.0.0-9cf.svg?style=flat)](https://snakemake.readthedocs.io) 

# HazGAN2

This repository contains a [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to generate spatially coherent climate multi-hazard event sets using extreme value theory and generative adversarial networks. The workflow is modular to facilitate new applications. It depends on a [Pytorch implementation]((https://github.com/NVlabs/stylegan2-ada-pytorch/issues/11)) of [StyleGAN2-ADA](https://arxiv.org/abs/2006.06676) with [differentiable augmentation](https://arxiv.org/abs/2006.10738) that provides stable training on small (~100 sample) datasets. 

The workflow is described in [this manuscript](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-3217/) and the rest of this README outlines basic usage.

## Branch management

The git branches are organised as follows:

- `main`: stable version
- `development`:
    - minor changes and bug fixes are added here first
    - feature branches are merged into here
    - periodically merge into `main` when stable
- `wip/<newfeature>`: branches for new features with potentially breaking changes


## Current status

### Checkpoints

- works with OUCE ERA5 reanalysis data ✅
- working with generalised Pareto marginal distributions ✅
- soft rescaling as a function of sample size (needs to be changed) ✅

### Pending changes

- rescaling as a function of reduced variate return levels
- using heavy-tailed latents to [match reduced variate tails better](https://arxiv.org/abs/2101.09113)
- implement [Weib–XIMIS](https://doi.org/10.3390/meteorology2030021) marginal fitting for more reliable wind speeds  
- switch to [Reiss and Thomas (2007)](https://doi.org/10.1007/978-3-7643-7399-3) tail dependence estimation

## Hardware and software

While the snakemake workflow can be run anywhere, StyleGAN2-ADA requires a CUDA-enabled NVIDIA GPU to train the GANs in a reasonable time. StyleGAN2-ADA has not been officially maintained, so it only works with CUDA versions up to 11.1. The recommended GPU is a NVIDIA V100. This workflow has been tested on the Oxford University Centre for the Environment's (OUCE) linux cluster with 1080Ti GPUs and on the University of Oxford's ARC cluster with V100 and RTX GPUs. See the official [Pytorch implementation]((https://github.com/NVlabs/stylegan2-ada-pytorch/issues/11)) for more information.

Because of the strict GPU requirements, it may be necessary to spread rules across different machines. The Snakefile and profiles have been hardcoded to work with OUCE and [ARC](https://www.arc.ox.ac.uk/getting-started-arc) computing clusters, but a new profile can be added for other machines, and the main Snakefile (in `workflow/Snakefile`) modified appropriately. See the `profiles/` directory for existing profiles.

Currently, the OUCE GPUs are not broken, so, while data processing is still done on the OUCE cluster, StyleGAN2 training is done on ARC.


## Quick start

To get started on a new machine, clone this repo and set up a micromamba environment with snakemake:

```bash
micromamba create -c conda-forge -c bioconda -n snakemake snakemake conda=24.7.1 -y

micromamba activate snakemake
# conda config --set-channel_priority strict  # (or change in .condarc file)

python -m pip install snakemake-executor-plugin-slurm # snakemake >= 9.0.0, if using SLURM

# set up the required conda environments
snakemake --profile profiles/cluster --conda-create-envs-only
snakemake --profile profiles/slurm --conda-create-envs-only

# dry run to view the workflow for current config
snakemake --profile profiles/local -n
```

## Things to do outside of this workflow

This workflow is only for generating event sets. To keep it tidy, downstream analysis should be done externally with code in a `projects/<project>/analysis` directory. This includes:

- Generating a `params.nc` file for the project
- All extra analysis (in `<project>/analysis`)

## Setup

### Apple silicon

For Apple Silicon, the R package `r-extremes` is not available on the conda `osx-arm64` subdirectory, so installation must be manually set to the `osx-64` subdirectory. If running a rule that will install the R environment for the first time, prefix the command with `CONDA_SUBDIR=osx-64`, e.g.,

```bash
CONDA_SUBDIR=osx-64 snakemake --profile profiles/local/ process_all_data --use-conda --cores 2
```

### OUCE cluster

When running for the first time, login nodes are extremely slow for creating conda environments. It's best to create the environments on an interactive compute node first:

```bash
srun -p Short --pty /bin/bash
micromamba activate snakemake
snakemake --profile profiles/cluster/ --conda-create-envs-only
```

After that, you need to run snakemake from the login node (as per current SoGE cluster configuration). To do this, it's best to open a screen session so that the job continues running if the connection drops:

```bash
screen -S snakemake
micromamba activate snakemake
cd path/to/hazGAN2
snakemake --profile profiles/slurm/ my_rule

# Ctrl+A then D to detach from screen session
# screen -r snakemake  # to reattach
```

How jobs are sent to SLURM is defined in the `config.yaml` file in each profile (local,cluster,arc,slurm). You can modify any of the profiles or make a new one to suit your device. To run the rule on an interactive compute node use the following command:

```bash
snakemake --profile profiles/cluster/ my_rule
```

or to send the job to the SLURM scheduler:

```bash
snakemake --profile profiles/slurm/ my_rule
```
### ARC cluster

The ARC cluster is the opposite to OUCE in that jobs must be submitted from a compute node, not the login node. To do this, use the `arc-submit.sh` script in the main directory. Modify this as suited. For now, ARC is only used for StyleGAN training.

## Further information

#### Basic modifications

To create a new project you need to make the following changes:

1. `config/config.yaml`: change `project` value
2. `config/projects/`: add a YAML file named `{myproject}.yaml` with the same structure as existing project YAMLs
3. `resources/params/`: use `resources/grids/era5.nc` to make `{myproject}.nc` in `resources/params/` with any spatial parameters for variable construction from raw ERA5 variables (see other param files for examples).

For each project, new Python and R functions can be added to the `src/` directory. See `projects/poweruk_winter` for an example.

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
