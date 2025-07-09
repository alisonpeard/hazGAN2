"""
Make parameters to use when pre-processing the ERA5 data.
"""
# %%
import os
from glob import glob
from tqdm import tqdm
import yaml
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def convert_360_to_180(ds):
    """Convert longitude from 0-360 to -180 to 180."""
    ds = ds.assign_coords(longitude=(ds.longitude + 180) % 360 - 180)
    return ds.sortby('longitude')

INPUT_DIR = "../../general_resources/params/"

# main
with open(os.path.join("..", "config.yaml"), "r") as stream:
    config = yaml.safe_load(stream)

dataset = config['dataset']
input_files = sorted(glob(os.path.join(INPUT_DIR, "*.nc")), reverse=True)
base = xr.open_dataset(os.path.join(INPUT_DIR, f"{dataset}.nc"))

xmin, xmax = config['longitude']['min'], config['longitude']['max']
ymin, ymax = config['latitude']['min'], config['latitude']['max']

base = convert_360_to_180(base)
params = base.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
params.sample.plot(ax=ax)
ax.add_feature(cfeature.COASTLINE)
ax.set_title("Sample Parameter")
ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
plt.show()

params.to_netcdf(os.path.join("..", "resources", "params.nc"),
                 encoding={
                     "sample": {
                         "dtype": "float32",
                         "zlib": True,
                         "complevel": 5
                     }
                 })

# %%
