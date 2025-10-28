# %%
import os
import sys
import yaml
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def format_poi(p):
    """Handle different formats of points of interest"""
    if isinstance(p, list):
        p = {"lon": p[0], "lat": p[1]}
    return p


def format_pois(p: list):
    """Format for extraction from xarray"""
    keys = list(p.keys())
    p = [format_poi(p_i) for p_i in p.values()]
    lats = [p_i["lat"] for p_i in p]
    lons = [p_i["lon"] for p_i in p]
    print(f"{lats=}, {lons=}")
    return keys, lats, lons


if __name__ == "__main__":
    # config and paths
    plt.rcParams["font.family"] = ["serif", "sans-serif", "monospace"][2]
    # bd = os.path.join("/hn01-home", "spet5107") # depends on device
    bd = os.path.join("/Users", "alison", "Local", "github")
    wd = "hazGAN2/projects/poweruk_winter"
    sys.path.append(os.path.join(bd, "hazGAN2", "workflow"))
    os.chdir(bd)

    with open(os.path.join(wd, "config.yaml"), "r") as stream:
        config = yaml.safe_load(stream)

    fields = list(config["fields"].keys())
    train_all = xr.open_dataset(os.path.join(wd, "results", "training", "data.nc"))
    print("Training data shapes:\n\n{}".format(train_all.sizes))

    pois, lats, lons = format_pois(config["points_of_interest"])
    num_pois = len(pois)

    lat_da = xr.DataArray(lats, dims="poi", coords={"poi": pois})
    lon_da = xr.DataArray(lons, dims="poi", coords={"poi": pois})
    train = train_all.sel(lat=lat_da, lon=lon_da, method="nearest")
    print("\nTrain for poi {}:\n\n{}".format(pois[0].title(), train.sel(poi=pois[0])))

    train.to_netcdf(os.path.join(wd, "results", "testing", "pois.nc"))
    print("Saved pois to file: {}".format(os.path.join(wd, "results", "testing", "pois.nc")))
    # %%
    