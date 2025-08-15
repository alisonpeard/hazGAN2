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



regions = ["East Midlands", "West Midlands",
           "South West England", "South Wales"]


with open(os.path.join("..", "config.yaml"), "r") as stream:
    config = yaml.safe_load(stream)

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

if config["event_subset"]["do"]:
    train["l2norm"] = np.sqrt(train.sel(field="u10").anomaly**2 + train.sel(field="v10").anomaly**2)
    mask = train.l2norm.max(dim=["lon", "lat"]) > 25
    train = train.isel(time=mask)
# %%
medians = train.medians.sel(month="December") # this should only have months for season

# %%
dataset = ["train", "generated"][0]
if dataset == "train":
    data = train.copy() #+ medians
else:
    data = gener.copy() #+ medians

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


# %% plot aggregate statistics
agg_func = ["mean", "median", "min", "max"][1]

data_agg = getattr(data, agg_func)(dim="time")
fig, axs = plt.subplots(1, 2, figsize=(10, 10), subplot_kw={"projection": ccrs.PlateCarree()})

ax = axs[0]
data_agg.vx.plot(
    ax=ax, cmap=cmo.speed,
    cbar_kwargs={"label": "mean gust speed (m/s)", "shrink": 0.6}
    )
data_agg.plot.streamplot(
    x='lon', y='lat', u='u10', v='v10',
    transform=ccrs.PlateCarree(), color="white", ax=ax, density=1.5,
    linewidth=0.5, arrowstyle='->', arrowsize=2
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.25)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_title("Mean wind speed and direction")

ax = axs[1]
data_agg.r30.plot(
    ax=ax, cmap=cmo.rain,
    cbar_kwargs={"label": "30-day rainfall (m)", "shrink": 0.6}
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.25)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_title("Mean 30-day antecedent rainfall")

fig.suptitle(f"{agg_func.capitalize()} {dataset.lower()}", fontsize=16, y=0.85)
plt.tight_layout()

# %% 
data_agg = data.std(dim="time")
fig, axs = plt.subplots(1, 3, figsize=(15, 10), subplot_kw={"projection": ccrs.PlateCarree()})

ax = axs[0]
data_agg.dx.plot(
    ax=ax, cmap=cmo.speed,
    vmin=60, vmax=120,
    cbar_kwargs={"label": "direction (m/s)", "shrink": 0.6}
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.25)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_title("Wind direction")

ax = axs[1]
data_agg.vx.plot(
    ax=ax, cmap=cmo.speed,
    vmin=2.5, vmax=5.5,
    cbar_kwargs={"label": "gust speed (m/s)", "shrink": 0.6}
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.25)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_title("Wind speed")

ax = axs[2]
data_agg.r30.plot(
    ax=ax, cmap=cmo.rain,
    vmin=0.02, vmax=0.14,
    cbar_kwargs={"label": "30-day rainfall (m)", "shrink": 0.6}
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.25)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_title("30-day antecedent rainfall")

fig.suptitle(f"Std {dataset.lower()}", fontsize=16, y=0.85)
plt.tight_layout()

# %% plot scatter plots for directon of specific locations in NE and SE
locations = {
    "North Sea": {"lon": 0.71512, "lat": 57.5237},
    "English Channel": {"lon": -0.2522, "lat": 50.3990}
}


fig, axs = plt.subplots(2, 2, figsize=(8, 7), sharex=True, sharey=True)

for i, (loc, coords) in enumerate(locations.items()):
    data_t = data.sel(lon=coords["lon"], lat=coords["lat"], method="nearest")
    ax = axs[0, i]
    speed = ax.scatter(data_t.u10, data_t.v10,
               c=data_t.vx, cmap=cmo.speed,
               s=10)
    fig.colorbar(speed, ax=ax, label="max gust speed (m/s)", shrink=0.6)
    
    ax = axs[1, i]
    dir = ax.scatter(data_t.u10, data_t.v10,
               c=data_t.dx, cmap="twilight_shifted",
               vmin=0, vmax=360,
               s=10)
    fig.colorbar(dir, ax=ax, label="direction (degrees)", shrink=0.6)
    
    for ax in axs[:, i]:
        ax.set_xlabel("u10")
        ax.set_ylabel("v10")
        ax.set_title(f"{loc}")
        ax.axhline(0, color='k', lw=0.5, alpha=0.1)
        ax.axvline(0, color='k', lw=0.5, alpha=0.1)

        m, b = np.polyfit(data_t.u10, data_t.v10, 1)
        ax.plot(data_t.u10, m * data_t.u10 + b,
                color="k", lw=0.75, alpha=0.75)
    
        ax.label_outer()


# %%
config["longitude"] # {'min': -10.875, 'max': 5.125}
config["latitude"] # {'min': 49.0, 'max': 64.75}

# %%
t = int(np.random.choice(range(data.sizes["time"]), 1)[0])

print(f"Selected time step: {t}")

data_t = data.isel(time=t)

resample = data_t.isel(lon=slice(None, None, 7), lat=slice(None, None, 7))

fig, axs = plt.subplots(1, 2, figsize=(18, 10), subplot_kw={"projection": ccrs.PlateCarree()})


ax = axs[0]
data_t.vx.plot(
    ax=ax, cmap=cmo.speed,
    cbar_kwargs={"label": "max gust speed (m/s)", "shrink": 0.6}
    )
data_t.plot.streamplot(
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
data_t.r30.plot(
    ax=ax, cmap=cmo.rain,
    cbar_kwargs={"label": "30-day rainfall (m)", "shrink": 0.6}
    )
gdf_subset.boundary.plot(color="k", ax=ax, linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.25)
ax.set_xlim(bbox[0], bbox[2])
ax.set_ylim(bbox[1], bbox[3])
ax.set_title("30-day antecedent rainfall")

fig.suptitle(f"{dataset.capitalize()} wind and rainfall for {t=}", fontsize=16, y=0.85)

plt.tight_layout()
# %%
"""
def calculate_nimgs(wildcards, years_of_samples=config["nyears"]):
    with zipfile.ZipFile(os.path.join(TRAINING_DIR, "images.zip"), 'r') as zip_ref:
        img_files = [f for f in zip_ref.namelist() if f.lower().endswith(('.png', '.jpg', '.jpeg', '.npy'))]
        nimgs = len(img_files)
    nyears = YEARN - YEAR0
    freq   = nimgs / nyears
    nsamples = int(freq * years_of_samples)
    print(f"Calculated number of samples: {nsamples} based on {nimgs} images over {nyears} years.")
    return nsamples

years_of_samples = 50
nyears = 2018 - 2004 + 1 = 15
nimgs = 429
freq = 429 / 15 = 28.6
nsamples = 28.6 * 50 = 1430
"""