"""
snakemake --profile local
snakemake --profile cluster
"""
import os
import yaml

configfile: "config.yaml"

env = config.get("environment", "local") 

YEARS = list(range(config["year0"], config["yearn"]))

INDIR = config['environments'][env]['era5dir']
WD = config['environments'][env]['datadir']
PROJDIR = os.path.join(WD, config['project'])
PROCESSING_DIR = os.path.join(PROJDIR, "processing")

CONDA = config['environments'][env]['environment']

FIELDS = config["fields"]
FIELD_LONG = config["field_names"]

for directory in [PROJDIR, PROCESSING_DIR]:
    os.makedirs(directory, exist_ok=True)

# load rules
include: "workflow/rules/get_data.smk"

