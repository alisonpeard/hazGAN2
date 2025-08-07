# %%
import os
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
from windrose import WindroseAxes


def direction(ds:xr.Dataset, u:str, v:str,
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


regions = ["East Midlands", "West Midlands",
           "South West England", "South Wales"]

path = os.path.join("..", "resources", "dno_license_areas.geojson")

gdf = gpd.read_file(path)
gdf_subset = gdf[gdf["region"].isin(regions)]

fig, ax = plt.subplots()
gdf.boundary.plot(color="k", ax=ax, linewidth=0.5)
gdf_subset.plot("region", ax=ax, categorical=True, legend=True)
bbox = gdf_subset.total_bounds
# %% Load real and generated events (and medians)
train = xr.open_dataset(os.path.join("..", "results", "training", "data.nc"))
gener = xr.open_dataset(os.path.join("..", "results", "generated", "netcdf", "data.nc"))

train["l2norm"] = np.sqrt(train.sel(field="u10").anomaly**2 + train.sel(field="v10").anomaly**2)
mask = train.l2norm.max(dim=["lon", "lat"]) > 25
train = train.isel(time=mask)
# %%
medians = train.medians.sel(month="December") # this should only have months for season

# %%
data = train.copy() + medians
# data = gener.copy() + medians

data["u10"] = data.sel(field="u10").anomaly
data["v10"] = data.sel(field="v10").anomaly
data["r30"] = data.sel(field="r30").anomaly

data = derive_variables(data)

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='windrose'))

# subset to bbox
data_wr = data.sel(
    lon=slice(bbox[0], bbox[2]),
    lat=slice(bbox[1], bbox[3])
    )

dirs = data_wr.dx.values.flatten()
speeds = data_wr.vx.values.flatten()

ax.bar(dirs, speeds, normed=True, opening=0.8, edgecolor='black',
        nsector=16, bins=10,
        cmap=cmo.speed, linewidth=0.5)
ax.set_legend()
# %%
t = np.random.choice(range(data.sizes["time"]), 1)[0]

print(f"Selected time step: {t}")

data = data.isel(time=t)

resample = data.isel(lon=slice(None, None, 7), lat=slice(None, None, 7))

fig, axs = plt.subplots(1, 2, figsize=(18, 10), subplot_kw={"projection": ccrs.PlateCarree()})


ax = axs[0]
data.vx.plot(
    ax=ax, cmap=cmo.speed,
    cbar_kwargs={"label": "max gust speed (m/s)", "shrink": 0.6}
    )
data.plot.streamplot(
    x='lon', y='lat', u='u10', v='v10', 
    transform=ccrs.PlateCarree(), color="white", ax=ax, density=1.5,
    linewidth=0.5, arrowstyle='->', arrowsize=2
    )

gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.5)

ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_xlim(bbox[0], bbox[2])
ax.set_ylim(bbox[1], bbox[3])
ax.set_title("Wind speed and direction")

ax = axs[1]
data.r30.plot(
    ax=ax, cmap=cmo.rain,
    cbar_kwargs={"label": "30-day rainfall (m)", "shrink": 0.6}
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_xlim(bbox[0], bbox[2])
ax.set_ylim(bbox[1], bbox[3])
ax.set_title("30-day antecedent rainfall")

plt.tight_layout()
# %%
