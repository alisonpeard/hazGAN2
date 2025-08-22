"""Plot histograms and return period plots."""
# %%
import os
import yaml
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


wd = os.path.join("..", "results")
# train_path = os.path.join(wd, "training", "data.nc")
train_path = os.path.join(wd, "generated", "netcdf", "train.nc")
gener_path = os.path.join(wd, "generated", "netcdf", "data.nc")

cfg_path = os.path.join("..", "config.yaml")

train = xr.open_dataset(train_path)
gener = xr.open_dataset(gener_path)

with open(cfg_path, "r") as stream:
    cfg = yaml.safe_load(stream)

# %%
def l2norm_max(ds:xr.Dataset, arg1:str, arg2:str,
              params={}, dim="time") -> xr.DataArray:
    l2norm = np.sqrt(ds[arg1]**2 + ds[arg2]**2)
    return l2norm.max(dim=dim)


if cfg["event_subset"]["do"]:
    print("Subsetting...")
    # func = cfg["event_subset"]["func"]
    func = l2norm_max
    args = cfg["event_subset"]["args"]
    value = cfg["event_subset"]["value"]

    for arg in args:
        train[arg] = train.sel(field=arg).anomaly
        gener[arg] = gener.sel(field=arg).anomaly   

    intensity_train = func(train, *args, dim=["lon", "lat"])
    intensity_gener = func(gener, *args, dim=["lon", "lat"])

    for arg in args:
        train = train.drop_vars(arg)
        gener = gener.drop_vars(arg)
    
    mask_train = intensity_train > value
    mask_gener = intensity_gener > value

    train = train.isel(time=mask_train)
    gener = gener.isel(time=mask_gener)
else:
    print("No subsetting applied.")

# %%
hist_kws = {"bins": 100, "density": True, "alpha": 0.5,
            "edgecolor": "k", "linewidth": 0.5}

for field, properties in cfg["fields"].items():
    title = properties.get("title", field)
    x = train.anomaly.sel(field=field).values.flatten()
    x_gen = gener.anomaly.sel(field=field).values.flatten()

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(x, label='Training', **hist_kws)
    ax.hist(x_gen, label='Generated', **hist_kws)
    ax.set_title(f"Histogram of {title.lower()} anomalies")
    ax.set_xlabel(f"{title} anomaly")
    ax.set_ylabel("Density")
    ax.legend()
    plt.show()

# %%
