"""Aggregate data for the DNO license areas of interest."""
# %%
import os
import yaml
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo


def direction(
        ds:xr.Dataset, u:str, v:str,
        params={}
        ) -> xr.DataArray:
    u = ds[u]
    v = ds[v]
    direction = (180 + 180 / np.pi * np.arctan2(v, u)) % 360
    return direction


def derive_variables(ds):
    ds["vx"] = np.sqrt(ds.u10**2 + ds.v10**2)
    ds["dx"] = direction(ds, "u10", "v10")
    return ds


regions = [
    "East Midlands", "West Midlands",
    "South West England", "South Wales"
    ]


with open(os.path.join("..", "config.yaml"), "r") as stream:
    config = yaml.safe_load(stream)

path = os.path.join("..", "resources", "dno_license_areas.geojson")

dno_regions_uk = gpd.read_file(path)
dno_regions = dno_regions_uk[dno_regions_uk["region"].isin(regions)]

fig, ax = plt.subplots()

dno_regions.boundary.plot(color="k", ax=ax, linewidth=0.5)
dno_regions.plot("region", ax=ax, categorical=True, legend=True)
bbox = dno_regions.total_bounds

# %% Load real and generated events (and medians)
train = xr.open_dataset(os.path.join("..", "results", "training", "data.nc"))
gener = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "data.nc"))

if config["event_subset"]["do"]:
    # subset to events with a max wind gust over a threshold
    train["l2norm"] = np.sqrt(train.sel(field="u10_gust").anomaly**2 + \
                              train.sel(field="v10_v10_gust").anomaly**2)
    mask = train.l2norm.max(dim=["lon", "lat"]) > config["event_subset"]["value"]
    train = train.isel(time=mask)

    gener = np.sqrt(gener.sel(field="u10_gust").anomaly**2 + \
                     gener.sel(field="v10_gust").anomaly**2)
    mask = gener.max(dim=["lon", "lat"]) > config["event_subset"]["value"]
    gener = gener.isel(time=mask)


# %% Work with the training or generated dataset
dataset = ["train", "generated"][1]

if dataset == "train":
    data = train.copy() #+ medians
else:
    data = gener.copy() #+ medians

data["u10"] = data.sel(field="u10_gust").anomaly
data["v10"] = data.sel(field="v10_gust").anomaly
data["r30"] = data.sel(field="r30").anomaly

data = derive_variables(data)
data["dx"] = data["dx"].__abs__()

# assign units to variables
data.vx.attrs["units"] = "m/s"
data.dx.attrs["units"] = "°"
data.r30.attrs["units"] = "mm"

# assign long names to variables
data.vx.attrs["long_name"] = "10m wind speed anomaly"
data.dx.attrs["long_name"] = "10m wind direction anomaly"
data.r30.attrs["long_name"] = "30-day accumulated precipitation anomaly"

# %% crop to license area
t = 100
agg_func = ["mean", "max"][1]

data_t = data.isel(time=t)

dno_region = dno_regions[dno_regions["region"] == "East Midlands"]
dno_bbox = dno_region.to_crs(27700).buffer(1500).to_crs(4326).total_bounds

lon_mask = (data_t.lon >= dno_bbox[0]) & (data_t.lon <= dno_bbox[2])
lat_mask = (data_t.lat >= dno_bbox[1]) & (data_t.lat <= dno_bbox[3])
data_region = data_t.isel(lon=lon_mask, lat=lat_mask)

fig, axs = plt.subplots(1, 3, figsize=(12, 4), subplot_kw={"projection": ccrs.PlateCarree()})

data_region.vx.plot(ax=axs[0], cmap=cmo.speed)
data_region.dx.plot(ax=axs[1], cmap=cmo.amp)
data_region.r30.plot(ax=axs[2], cmap=cmo.rain)

for ax in axs:
    dno_region.boundary.plot(color="k", ax=ax, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

# %% aggregate over license area
from regionmask import mask_geopandas

lon = data_region.lon.values
lat = data_region.lat.values
dno_mask = mask_geopandas(dno_region, lon, lat)
# %% visualise mask
data_region["dno_mask"] = dno_mask

fig, axs = plt.subplots(1, 3, figsize=(12, 4), subplot_kw={"projection": ccrs.PlateCarree()})

data_region.vx.plot(ax=axs[0], cmap=cmo.speed)
data_region.dx.plot(ax=axs[1], cmap=cmo.amp)
data_region.r30.plot(ax=axs[2], cmap=cmo.rain)

for ax in axs:
    dno_region.boundary.plot(color="k", ax=ax, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

    data_region.dno_mask.plot(alpha=0.5, ax=ax)

# %% apply to all data
data["dno_mask"] = mask_geopandas(dno_regions[dno_regions["region"] == "East Midlands"], data.lon, data.lat)
data_agg = data.groupby("dno_mask").reduce(getattr(np, agg_func))

# %%
fig, ax = plt.subplots(1, 3, figsize=(12, 4))
data_agg.vx.plot.hist(ax=ax[0], bins=20)
data_agg.dx.plot.hist(ax=ax[1], bins=20)
data_agg.r30.plot.hist(ax=ax[2], bins=20)
# %%
data_agg
# %% load the coefs
import pandas as pd
coefs = pd.read_csv(os.path.join("..", "resources/E_Mid_All_Variable_Model_coefficients_wind-thresh=20.csv"))
coefs
# %%
# data_agg must be parsed to items [1, vx, r30, dx, spring, summer, autumn] with last three all zero

X = np.hstack([
    np.ones(data_agg.vx.shape).reshape(-1, 1),
    data_agg.vx.values.reshape(-1, 1),
    data_agg.r30.values.reshape(-1, 1),
    data_agg.dx.values.reshape(-1, 1),
    np.zeros(data_agg.vx.shape).reshape(-1, 1),
    np.zeros(data_agg.vx.shape).reshape(-1, 1),
    np.zeros(data_agg.vx.shape).reshape(-1, 1)
    ])

X.shape

# %%
# logistic regression model using coefs
from scipy.special import expit

y = expit(np.dot(X, coefs["Coefficients"].values))
# %%
plt.hist(y, bins=20)
# %%
