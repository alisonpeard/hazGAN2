"""
snakemake --profile local
snakemake --profile cluster
"""
import os
import yaml

configfile: "config.yaml"

env = config.get("environment", "local") 

project_config = "projects/{}.yaml".format(config["project"])
configfile: project_config

YEARS = list(range(config["year0"], config["yearn"]))

INDIR = config['environments'][env]['era5dir']
WD = config['environments'][env]['datadir']
PROJDIR = os.path.join(WD, "hazGAN", config['project'])
PROCESSING_DIR = os.path.join(PROJDIR, "processing")

CONDA = config['environments'][env]['environment']

FIELDS = config["fields"]

for directory in [PROJDIR, PROCESSING_DIR]:
    os.makedirs(directory, exist_ok=True)

# load rules
include: "workflow/rules/get_data.smk"

