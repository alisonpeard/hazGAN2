 # hazGAN2 Readme

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
│       │   ├── input/
│       │   └── generated/
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