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

from mangroves.model import mangroveDamageModel


MONTH = 9

if __name__ == "__main__":
    # load project configuration
    with open(os.path.join("..", "config.yaml"), 'r') as stream:
        config = yaml.safe_load(stream)

    # load generated data
    THRESHOLD = config['event_subset']["value"]
    NYRS = config["yearn"] - config["year0"]

    train = xr.open_dataset(os.path.join("..", "results", "training", "data.nc"))
    gener = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "data.nc"))
    indep = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "independent.nc"))
    depen = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "dependent.nc"))
    medians = train["medians"]

    #  extract the extreme samples
    if config["event_subset"] is not None:
        subset = config["event_subset"]
        threshold = config["event_subset"]['value']
        train['intensity'] = getattr(train.sel(field=subset["field"]).anomaly, subset["func"])(dim=['lon', 'lat'])
        mask = (train['intensity'] > subset["value"]).values
        idx  = np.where(mask)[0]
    else:
        idx = np.arange(train.sizes["time"])

    nevents = train.sizes["time"]
    nextremes = len(idx)

    train = train.isel(time=idx)
    medians = medians.isel(month=MONTH)

    # get event sets
    train["value"] = train["anomaly"] + medians
    gener["value"] = gener["anomaly"] + medians
    indep["value"] = indep['anomaly'] + medians
    depen["value"] = depen['anomaly'] + medians

    # check for nans
    for ds in [train, gener, indep, depen]:
        assert ds.isnull().sum() == 0, "Nans found in dataset"

    # predict mangrove damages
    model = mangroveDamageModel()
    print(model)

    train_damages = model.predict(train, ["value"])
    gener_damages = model.predict(gener, ["value"])
    indep_damages = model.predict(indep, ["value"])
    depen_damages = model.predict(depen, ["value"], first_dim="rp")

    # rename the <>_damages to damage_prob
    rename_dict = {"value_damage": "damage_prob",}
    train_damages = train_damages.rename(rename_dict)
    gener_damages = gener_damages.rename(rename_dict)
    indep_damages = indep_damages.rename(rename_dict)
    depen_damages = depen_damages.rename(rename_dict)

    # add rate attribut to each dataset
    if config["event_subset"] is not None:
        train_damages.attrs["rate"] = nextremes / NYRS
    else:
        train_damages.attrs["rate"] = nevents / NYRS
    gener_damages.attrs["rate"] = nextremes / NYRS
    indep_damages.attrs["rate"] = nevents / NYRS
    # depen_damages.attrs["rate"] = None

    # make return period maps from fake winds
    outdir = os.path.join("..", "results", "mangroves", "damage_fields")
    os.makedirs(outdir, exist_ok=True)

    train_damages.to_netcdf(os.path.join(outdir, "train.nc"))
    gener_damages.to_netcdf(os.path.join(outdir, "gener.nc"))
    indep_damages.to_netcdf(os.path.join(outdir, "indep.nc"))
    depen_damages.to_netcdf(os.path.join(outdir, "depen.nc"))

    # plot somthin'
    fig, axs = plt.subplots(1, 4, figsize=(16, 4))
    train_damages["damage_prob"].plot(ax=axs[0])
    gener_damages["damage_prob"].plot(ax=axs[1])
    indep_damages["damage_prob"].plot(ax=axs[2])
    depen_damages["damage_prob"].plot(ax=axs[3])

# %%