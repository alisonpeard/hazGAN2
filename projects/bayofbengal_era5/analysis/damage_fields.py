"""
Step 1: Mangrove damage probability fields
convert ERA5 data, generated samples and benchmark data to netcdf
load samples as netcdf
add monthly medians
load traning data -> 
get yearly rate
define damage ufunc and apply it to each netcdf samples, train, test
NEW: save damages
"""
# %%
import os
import yaml
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .mangroves.model import mangroveDamageModel
# from analysis import load_samples, get_monthly_medians, λ

if __name__ == "__main__":
    # load project configuration
    with open(os.path.join("..", "config.yaml"), 'r') as stream:
        config = yaml.safe_load(stream)

    # load generated data
    THRESHOLD = config['event_subset']["value"]
    NYRS = config["yearn"] - config["year0"]
    MONTH = 9

    train = xr.open_dataset(os.path.join("..", "results", "training", "data.nc"))
    gener = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "data.nc"))
    indep = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "independent.nc"))
    depen = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "dependent.nc"))
    
    # load and process monthly medians (deseaonalize)
    medians = pd.read_csv(os.path.join("..", "results", "processing", "medians.csv"))
    medians = medians[medians["month"] == MONTH]
    medians = medians.set_index(["lat", "lon"])
    medians = medians.to_xarray()

    # extract the extreme samples
    if config["event_subset"] is not None:
        subset = config["event_subset"]
        threshold = config["event_subset"]['threshold']
        train['intensity'] = getattr(train.sel(field=subset["field"]).anomaly, subset["func"])(dim=['lon', 'lat'])
        mask = (train['intensity'] > subset["value"]).values
        idx  = np.where(mask)[0]
    else:
        idx = np.arange(train.sizes["time"])

    nevents = train.sizes["time"]
    nextremes = len(idx)

    # get event sets
    train["value"] = train["anomaly"] + medians
    gener["value"] = gener["anomaly"] + medians
    indep["value"] = indep['anomaly'] + medians
    depen["value"] = depen['anomaly'] + medians

    # convert dependent uniform array to a 1-d array
    assert np.isclose(depen['uniform'].std(dim=['lat', 'lon']).max(), 0)
    depen['uniform'] = depen['uniform'].mean(dim=['lat', 'lon'])

    # check for nans
    for ds in [train, gener, indep, depen]:
        assert ds.isnull().sum() == 0, "Nans found in dataset"
        # NaNs found here! Why?

    # predict mangrove damages
    model = mangroveDamageModel()
    train_damages = model.predict(train, ["train"])
    gener_damages = model.predict(gener, ["gener"])
    indep_damages = model.predict(indep, ["indep"])
    depen_damages = model.predict(depen, ["depen"])

    # rename the <>_damages to damage_prob
    train_damages = train_damages.rename({"train_damage": "damage_prob"})
    gener_damages = gener_damages.rename({"gener_damage": "damage_prob"})
    indep_damages = indep_damages.rename({"indep_damage": "damage_prob"})
    depen_damages = depen_damages.rename({"depen_damage": "damage_prob"})

    # add rate attribut to each dataset
    train_damages.attrs["rate"] = nevents / NYRS
    gener_damages.attrs["rate"] = nextremes / NYRS
    indep_damages.attrs["rate"] = nevents / NYRS
    depen_damages.attrs["rate"] = None

    # make return period maps from fake winds
    outdir = os.makedirs(os.path.join("..", "results", "mangroves"), exist_ok=True)
    outdir = os.makedirs(os.path.join(outdir, "damage_fields"), exist_ok=True)

    train_damages.to_netcdf(os.path.join(outdir, "train.nc"))
    gener_damages.to_netcdf(os.path.join(outdir, "gener.nc"))
    indep_damages.to_netcdf(os.path.join(outdir, "indep.nc"))
    depen_damages.to_netcdf(os.path.join(outdir, "depen.nc"))


# %%