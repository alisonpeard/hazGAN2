"""
Forum: https://forum.ecmwf.int/t/unusual-10m-wind-distribution-over-bay-of-bengal/7399 

Solutions:

"""
# %%
import os
from glob import glob
from tqdm import tqdm
import yaml
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

INPUT_DIR = "/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/"
PROBLEM_DATES = [
    "1952-05-09",
    "1992-04-15",
    "1997-05-16"
]


# main
with open(os.path.join("..", "config.yaml"), "r") as stream:
    config = yaml.safe_load(stream)

xmin, xmax = config['longitude']['min'], config['longitude']['max']
ymin, ymax = config['latitude']['min'], config['latitude']['max']

input_files = sorted(glob(os.path.join(INPUT_DIR, "*.nc")), reverse=True)

base = xr.open_dataset(os.path.join("..", "..", "general_resources", "params", "era5.nc"))
params = base.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))
# params.sample.plot()
params.to_netcdf(os.path.join("..", "resources", "params", "era5_bayofbengal.nc"),
                 encoding={
                     "sample": {
                         "dtype": "float32",
                         "zlib": True,
                         "complevel": 5
                     }
                 })
