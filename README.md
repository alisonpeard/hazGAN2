# hazGAN2
Centralised repo for hazGAN methods and snakemake workflows, will publish later.

Repository structure:
```
hazGAN2/
├── .env
├── README.md
├── config
│   ├── csv files with metadata e.g., cost curves
│   ├── folder with arc cluster config
│   └── config.yaml: paths to datasets
├── StyleGAN2-DA
│   ├── environments
│   ├── pyproject.toml
│   └── src
│       └── *
├── hazGAN
│   ├── environments
│   ├── pyproject.toml
│   └── src
│       └── *
└── scripts
    ├── data_acquisition
    ├── data_processing
    ├── training
    ├── validation
    └── mangroves
```

Data output structure:
```
data/
├── bayofbengal
│   ├── processing
│   ├── training
│   │   └── 64x64
│   │   │   ├── input
│   │   │   └── generated
│   └── mangroves
└── unitedkingdom
    ├── processing
    ├── training
    └── analysis
 ```

# ERA5 Climate Analysis Project Structure

## Directory Structure

```
era5-project/
├── Snakefile                # Main Snakefile that orchestrates the workflow
├── config.yaml              # Default configuration
├── config/                  # Machine-specific configurations
│   ├── machine1.yaml
│   ├── machine2.yaml
│   ├── hpc.yaml
│   └── slurm/               # SLURM profile for cluster execution
│       └── config.json
├── workflow/
│   ├── rules/               # Modular rule files
│   │   ├── preprocessing.smk
│   │   ├── analysis.smk
│   │   └── visualization.smk
│   └── scripts/             # Python scripts called by rules
│       ├── process_era5.py
│       ├── combine_years.py
│       ├── extract_variable.py
│       ├── create_climatology.py
│       ├── trend_analysis.py
│       ├── extreme_events.py
│       ├── seasonal_patterns.py
│       ├── spatial_correlation.py
│       ├── teleconnection_analysis.py
│       ├── plot_trend_map.py
│       ├── plot_extreme_events.py
│       ├── plot_seasonal_cycle.py
│       ├── plot_timeseries.py
│       ├── plot_correlation_heatmap.py
│       ├── create_dashboard.py
│       └── config_utils.py  # Helper for config inheritance
├── logs/                    # Log files
├── sbatch_dump/             # SLURM output files
└── profiles/                # Additional execution profiles
    └── standard/            # Standard configuration
        └── config.yaml
```

## Snakemake Best Practices

### 1. Use a Modular Structure

- Split your workflow into logical modules (preprocessing, analysis, visualization)
- Use `include:` directives to bring in rule files
- Keep the main Snakefile simple and focused on workflow organization

### 2. Define Clear Target Rules

Instead of a single `all` rule, define multiple target rules that represent specific pipeline stages:

- `rule preprocess:` - Data preparation only
- `rule analyze:` - Run analyses on preprocessed data
- `rule visualize:` - Generate visualizations
- `rule full_pipeline:` - Complete end-to-end workflow

### 3. Use Wildcards Effectively

- Wildcards allow parameter-based execution (e.g., `{year}`, `{variable}`)
- Define `wildcard_constraints:` to restrict wildcard values
- Use `expand()` to generate combinations of wildcard values

### 4. Document Your Rules

- Add docstrings to every rule explaining its purpose
- Include parameters, input/output descriptions
- Add comments for complex logic

### 5. Manage Resources Appropriately

- Specify resource requirements for each rule
- Adjust based on task complexity (more memory for larger datasets)
- Use different resource specifications for different execution environments

### 6. Implement Error Handling

- Use the `onsuccess:` and `onerror:` directives for notifications
- Comprehensive logging in each rule
- Consider retry mechanisms for unreliable operations

### 7. Version Control

- Keep your Snakemake workflow in version control (git)
- Include `.gitignore` for generated files
- Document changes in a changelog

### 8. Testing

- Create a small test dataset
- Define test rules that run quickly
- Verify outputs against expected results

## Running Specific Pipeline Parts

### Run specific years

```bash
snakemake --config year_start=1980 year_end=1990
```

### Run specific analyses

```bash
snakemake trend_analysis extreme_events
```

### Generate specific visualizations

```bash
snakemake --forcerun trend_maps seasonal_cycle
```

### Rerun specific steps for specific years

```bash
snakemake --touch preprocess
snakemake process_era5_year --config years=[1980,1981,1982]
```

### Run with increased parallelism

```bash
snakemake --profile config/slurm --jobs 100 full_pipeline
```

## Advanced Features

### Checkpoints

For dynamic dependencies (when outputs are determined during execution):

```python
checkpoint process_dataset:
    # checkpoint definition
    
rule analyze_results:
    input:
        lambda wildcards: get_checkpoint_outputs(wildcards)
```

### DAG Visualization

Visualize your workflow as a directed acyclic graph:

```bash
snakemake --dag | dot -Tsvg > workflow.svg
```

### Report Generation

Generate a comprehensive HTML report:

```bash
snakemake --report report.html
```

### Remote File Support

Access files from remote storage (S3, Google Cloud, etc.):

```python
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

rule download_data:
    input:
        S3.remote("bucket/path/to/data.nc")
    output:
        "local/path/data.nc"
```

## Integration with Other Tools

- **Conda**: Use conda environments for each rule
- **Singularity/Docker**: Package dependencies in containers
- **GitHub Actions**: CI/CD for workflow testing
- **DVC**: Data version control alongside workflow



## Recommended directory structure from docs
```
├── .gitignore
├── README.md
├── LICENSE.md
├── workflow
│   ├── rules
|   │   ├── module1.smk
|   │   └── module2.smk
│   ├── envs
|   │   ├── tool1.yaml
|   │   └── tool2.yaml
│   ├── scripts
|   │   ├── script1.py
|   │   └── script2.R
│   ├── notebooks
|   │   ├── notebook1.py.ipynb
|   │   └── notebook2.r.ipynb
│   ├── report
|   │   ├── plot1.rst
|   │   └── plot2.rst
|   └── Snakefile
├── config
│   ├── config.yaml
│   └── some-sheet.tsv
├── results
└── resources
```