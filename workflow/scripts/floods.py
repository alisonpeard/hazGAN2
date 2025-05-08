# %%
import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from cartopy import crs as ccrs


ds = xr.open_dataset("/Users/alison/Documents/DPhil/github.nosync/hazGAN2/results/caribbeanrisk/training/data.nc")
# ds = xr.open_dataset("/Users/alison/Documents/DPhil/github.nosync/hazGAN2/results/caribbeanrisk/processing/data_all.nc")
ds
ds['sro'] = ds['anomaly'].sel(field="sro")
ds['tp'] = ds['anomaly'].sel(field="tp")
ds['ws'] = ds['anomaly'].sel(field="ws")

ds["maxsro"] = ds['sro'].sum(dim=["lon", "lat"])
ds = ds.sortby("maxsro", ascending=False)
# %%
fig, axs = plt.subplots(1, 3, figsize=(16, 4))

ds['tp'].isel(time=0).plot(cmap='Blues', ax=axs[0])
ds['sro'].isel(time=0).plot(cmap='Blues', ax=axs[1])
ds['ws'].isel(time=0).plot(cmap='Spectral_r', ax=axs[2])
# %%
from shapely.geometry import box
import geopandas as gpd

datadir = "/Users/alison/Documents/DPhil/data/hydrobasins/hybas_as_lev01-12_v1c"

LEVEL = 6
# 10 is too grainy-looking
# 5 looks nice but pretty lowres
# 8 has those gaps..
EXTENT = box(80, 10, 95, 25)
file = f"hybas_as_lev{str(LEVEL).zfill(2)}_v1c.shp"
gdf = gpd.read_file(os.path.join(datadir, file))
gdf = gdf.to_crs("EPSG:4326")
basins = gdf.clip(EXTENT)

# %%
dem = xr.open_dataset("/Users/alison/Documents/DPhil/data/etopo1/etopo1-bay-of-bengal.nc")

def accumulate_to_basins(sro:xr.DataArray, basins:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Calculate inundation from sro and basins
    """
    # copy incase they are mutable
    sro = sro.copy()
    basins = basins.copy()

    sro = sro.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
    sro = sro.rio.write_crs("epsg:4326")
    sro = sro.rio.clip(basins.geometry, basins.crs, drop=True)
    sro = sro.to_dataframe().reset_index()[['time', 'lon', 'lat', 'sro']]
    sro = gpd.GeoDataFrame(sro, geometry=gpd.points_from_xy(sro.lon, sro.lat))
    sro = sro.set_crs("EPSG:4326")
    sro = sro.dropna(subset='sro')

    sro = gpd.sjoin(basins, sro, how="left", predicate="intersects")
    sro = sro.drop(columns=['index_right', 'lon', 'lat'])
    sro = sro[["HYBAS_ID","NEXT_DOWN", "SUB_AREA", "SORT", "geometry", "sro", "time"]]
    sro = sro.groupby(["HYBAS_ID","NEXT_DOWN", "SUB_AREA", "SORT", "geometry", "time"]).mean()
    sro = sro.reset_index()

    # assumes relatively homogeneous sro across basin
    sro = gpd.GeoDataFrame(sro, geometry=sro.geometry)
    sro["volume"] = sro["sro"] * sro["SUB_AREA"]

    return sro

sro_basins = accumulate_to_basins(ds['sro'], basins)

t = sro_basins["time"].unique().to_numpy()

fig, ax = plt.subplots(2, 5, figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree()})
for i in range(5):
    sro_today = sro_basins[sro_basins["time"] == t[i]]
    sro_today.plot("volume", cmap="Blues", edgecolor='k',
                   linewidth=0.1, ax=ax[0, i],
                   legend=True,
                   legend_kwds={'label': "Runoff (m3)", 'orientation': "horizontal"})
    ds.sel(time=t[i]).tp.plot.contourf(ax=ax[1, i], levels=10, cmap='Blues', alpha=.4, linewidths=.5,
                                        add_colorbar=True, vmin=0,
                                        cbar_kwargs={'label': "Precipitation (m)", 'orientation': "horizontal"})

# %%
def propagate_runoff(sro_basins:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Propagate runoff upstream to downstream basins
    """
    # 1 is most downstream
    sro_basins = sro_basins.set_index("HYBAS_ID")
    flow_accumulator = sro_basins.copy()

    i = 0
    while (len(flow_accumulator) > 1):
        # print(i, len(flow_accumulator))
        flow_accumulator = flow_accumulator[["NEXT_DOWN", "volume"]].groupby(["NEXT_DOWN"])
        flow_accumulator = flow_accumulator.agg(inflow=("volume", "sum")).reset_index()
        flow_accumulator = flow_accumulator.rename(columns={"NEXT_DOWN": "HYBAS_ID"})
        flow_accumulator = flow_accumulator.set_index(["HYBAS_ID"])

        # add new inflows to basin volumes
        sro_basins = sro_basins.merge(flow_accumulator[["inflow"]], how="left", left_index=True, right_index=True)
        sro_basins["inflow"] = sro_basins["inflow"].fillna(0)
        sro_basins["volume"] = sro_basins["volume"] + sro_basins["inflow"]
        sro_basins = sro_basins.drop(columns=["inflow"])

        # assign updated basin volumes to flow accumulator
        flow_accumulator = flow_accumulator.merge(sro_basins[["NEXT_DOWN", "volume"]], how="left", left_index=True, right_index=True)
        i += 1

    return sro_basins

# sro_basins_day0 = sro_basins[sro_basins["time"] == sro_basins["time"].min()]

runoffs = []
for i in range(len(t)):
    print(f"Processing time {i} out of {len(t)}")
    sro_today = sro_basins[sro_basins["time"] == t[i]]
    accumulated_runoff = propagate_runoff(sro_today)
    runoffs.append(accumulated_runoff)

runoffs = pd.concat(runoffs)
runoffs = gpd.GeoDataFrame(runoffs, geometry=runoffs.geometry)

# %%
vmin = sro_basins["volume"].min()
vmax = sro_basins["volume"].max()

fig, ax = plt.subplots(4, 5, figsize=(16, 12), subplot_kw={'projection': ccrs.PlateCarree()})
for i in range(5):
    sro_today = sro_basins[sro_basins["time"] == t[i]]
    runoff_today = runoffs[runoffs["time"] == t[i]]
    wind_today = ds.sel(time=t[i]).ws
    tp_today = ds.sel(time=t[i]).tp

    wind_today.plot(ax=ax[0, i], levels=10, cmap='Spectral_r', linewidths=.5,
                                        add_colorbar=True, vmin=0,
                                        cbar_kwargs={'label': "Wind speed (mps)", 'orientation': "horizontal"}
    )

    sro_today.plot("volume", cmap="Blues", edgecolor='k',
                   linewidth=0.1, ax=ax[2, i],
                   legend=True, vmin=vmin, vmax=vmax,
                   legend_kwds={'label': "Runoff (m3)", 'orientation': "horizontal"})


    runoff_today.plot("volume", cmap="Blues", edgecolor='k',
                   linewidth=0.1, ax=ax[3, i],
                   legend=True, vmin=vmin, vmax=vmax,
                   legend_kwds={'label': "Runoff (m3)", 'orientation': "horizontal"})
    
    tp_today.plot.contourf(ax=ax[1, i], levels=10, cmap='Blues', linewidths=.5,
                                        add_colorbar=True, vmin=0,
                                        cbar_kwargs={'label': "Precipitation (m)", 'orientation': "horizontal"})
plt.tight_layout()
# %%
fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={'projection': ccrs.PlateCarree()})
dem.elevation.plot(cmap='viridis', vmin=0, add_colorbar=False, ax=ax)
# %% OLD (for reference)

# import rioxarray

# ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
# ds.rio.write_crs("epsg:4326", inplace=True)
# clipped = ds.rio.clip(gdf.geometry, gdf.crs, drop=True)
# runoff = clipped['sro'].to_dataframe().reset_index()[['time', 'lon', 'lat', 'sro']]
# runoff = gpd.GeoDataFrame(runoff, geometry=gpd.points_from_xy(runoff.lon, runoff.lat))
# runoff = runoff.set_crs("EPSG:4326")
# runoff = runoff.dropna(subset='sro')

# for i in range():
# fig, ax = plt.subplots(1, 5, figsize=(10, 10))
# runoff.plot(ax=ax, column='sro', cmap='Blues', legend=True)
# basins.plot(ax=ax, color='none', edgecolor='black')

# intersection = gpd.sjoin(basins, runoff, how="left", predicate="contains")
# intersection = intersection.drop(columns=['index_right', 'lon', 'lat'])
# intersection = intersection[["HYBAS_ID","NEXT_DOWN", "SUB_AREA", "SORT", "geometry", "sro"]].groupby(["HYBAS_ID","NEXT_DOWN", "SUB_AREA", "SORT", "geometry"]).mean()
# intersection = intersection.reset_index()

# intersection = gpd.GeoDataFrame(intersection, geometry=intersection.geometry)

# 1 is most downstream
# intersection = intersection.sort_values(by="SORT", ascending=False)
# intersection["volume"] = intersection["sro"] * intersection["SUB_AREA"]
# intersection.plot("volume", cmap="Blues", legend=True)


# from matplotlib.colors import LogNorm
# VOI = "volume" # "volume" or "sro"
# fig, axs = plt.subplots(1, 3, figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree()})

# # propagate sro upstream --> downstream by summation
# totals = intersection.copy().set_index("HYBAS_ID")
# totals.sort_values(by=VOI, ascending=False).head()
# totals.plot(VOI, cmap="Blues", legend=True, ax=axs[0])
# flow_accumulator = totals.copy()

# i = 0
# while (len(flow_accumulator) > 1):
# # while i < 3:
#     print(i, len(flow_accumulator))
#     flow_accumulator = flow_accumulator[["NEXT_DOWN", VOI]].groupby("NEXT_DOWN").agg(inflow=(VOI, "sum")).reset_index()
#     flow_accumulator = flow_accumulator.rename(columns={"NEXT_DOWN": "HYBAS_ID"})
#     flow_accumulator = flow_accumulator.set_index("HYBAS_ID")

#     # add new inflows to basin volumes
#     totals = totals.merge(flow_accumulator[["inflow"]], how="left", left_index=True, right_index=True)
#     totals["inflow"] = totals["inflow"].fillna(0)
#     totals[VOI] = totals[VOI] + totals["inflow"]
#     totals = totals.drop(columns=["inflow"])

#     # assign updated basin volumes to flow accumulatoe
#     flow_accumulator = flow_accumulator.merge(totals[["NEXT_DOWN", VOI]], how="left", left_index=True, right_index=True)
#     i += 1

# totals.plot(VOI, cmap="Blues", norm=LogNorm(), legend=True, ax=axs[1]) # 
# dem.elevation.plot(cmap='viridis', vmin=0, ax=axs[2], add_colorbar=False)
# # same plot but log colorscale
# totals.sort_values(by=VOI, ascending=False).head()
# %%
