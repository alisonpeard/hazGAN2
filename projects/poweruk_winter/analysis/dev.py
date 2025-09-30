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
from windrose import WindroseAxes #noqa: E402


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
    train["l2norm"] = np.sqrt(train.sel(field="u10_gust").anomaly**2 + \
                              train.sel(field="v10_v10_gust").anomaly**2)
    mask = train.l2norm.max(dim=["lon", "lat"]) > config["event_subset"]["value"]
    train = train.isel(time=mask)

    gener = np.sqrt(gener.sel(field="u10_gust").anomaly**2 + \
                     gener.sel(field="v10_gust").anomaly**2)
    mask = gener.max(dim=["lon", "lat"]) > config["event_subset"]["value"]
    gener = gener.isel(time=mask)

# %%
import pandas as pd

fields = ["u10_gust", "v10_gust", "r30"]

medians_file = os.path.join("..", "results", "processing", "medians.parquet")
medians = pd.read_parquet(medians_file)
medians = medians.astype({"lat": float, "lon": float})

medians = medians.sort_values(["month", "lat", "lon"], ascending=[True, False, True])
print(medians.head())
medians = medians[fields].values.reshape([12, 64, 64, 3])
plt.imshow(medians[:,:,0,0], cmap=cmo.speed);


# %%
train.sel(field="u10_gust").medians.isel(month=0).plot(cmap=cmo.speed)
# %%
medians = train.medians.sel(month="December") # this should only have months for season
medians.sel(field="u10_gust").plot(cmap=cmo.speed)
# %%
medians = medians.to_dataset(name="medians")
medians["u10"] = medians.sel(field="u10_gust").medians
medians["v10"] = medians.sel(field="v10_gust").medians
medians["r30"] = medians.sel(field="r30").medians
medians = derive_variables(medians)

# %%
medians.dx.plot(cmap=cmo.balance)
# %% Work with the training or generated dataset
dataset = ["train", "generated"][0]

if dataset == "train":
    data = train.copy() #+ medians
else:
    data = gener.copy() #+ medians

data["u10"] = data.sel(field="u10_gust").anomaly
data["v10"] = data.sel(field="v10_gust").anomaly
data["r30"] = data.sel(field="r30").anomaly

data = derive_variables(data)

# %% crop to license area
t = 0
agg_func = ["mean", "max"][1]

data_t = data.isel(time=t)

dno_region = dno_regions[dno_regions["region"] == "East Midlands"]
dno_bbox = dno_region.to_crs(27700).buffer(1500).to_crs(4326).total_bounds

lon_mask = (data_t.lon >= dno_bbox[0]) & (data_t.lon <= dno_bbox[2])
lat_mask = (data_t.lat >= dno_bbox[1]) & (data_t.lat <= dno_bbox[3])
data_region = data_t.isel(lon=lon_mask, lat=lat_mask)

fig, axs = plt.subplots(1, 3, figsize=(8,2), subplot_kw={"projection": ccrs.PlateCarree()})

ax = axs[0]
data_region.vx.plot(ax=axs[0], cmap=cmo.speed)
data_region.dx.plot(ax=axs[1], cmap=cmo.amp)
data_region.r30.plot(ax=axs[2], cmap=cmo.rain)

for ax in axs:
    dno_region.boundary.plot(color="k", ax=ax, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.set_title(f"{agg_func.capitalize()} wind speed over {dno_region.iloc[0]['region']}")



#