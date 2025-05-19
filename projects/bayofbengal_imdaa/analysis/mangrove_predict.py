# Step 1: Mangrove damage probability fields
# convert samples to netcdf and add dependence assumption fields
# load samples as netcdf (including dependentce assumptions)
# add monthly medians
# load traning data -> 
# get yearly rate
# define damage ufunc and apply it to each netcdf samples, train, test
# NEW: save damages

# Step 2: Intersect damage fields with mangrove fields
# load mangroves
# intersect mangroves with damage fields
# get mangrove damages
# Calculate damagearea

# %%
import os
import yaml
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from hazGAN.mangrove_demo import mangroveDamageModel
# from analysis import load_samples, get_monthly_medians, λ

if __name__ == "__main__":
    # load project configuration
    with open(os.path.join("..", "config.yaml"), 'r') as stream:
        config = yaml.safe_load(stream)

    # load generated data
    THRESHOLD = config['event_subset']['threshold']

    train = xr.open_dataset(os.path.join("..", "results", "training", "data.nc"))
    gener = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "data.nc"))
    indep = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "independent.nc"))
    depen = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "dependent.nc"))
    # medians = get_monthly_medians(data_dir, "September") # should be size 1 x 64 x 64 x 3

    nobs = train.sizes["time"]
    nyrs = config["yearn"] - config["year0"]

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

    # make return period maps from fake winds
    outdir = os.makedirs(os.path.join("..", "results", "mangroves"), exist_ok=True)
    train_damages.to_netcdf(os.path.join(outdir, "train_damages.nc"))
    gener_damages.to_netcdf(os.path.join(outdir, "gener_damages.nc"))
    indep_damages.to_netcdf(os.path.join(outdir, "indep_damages.nc"))
    depen_damages.to_netcdf(os.path.join(outdir, "depen_damages.nc"))

# %%