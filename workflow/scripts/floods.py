# %%
import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt


ds = xr.open_dataset("/Users/alison/Documents/DPhil/github.nosync/hazGAN2/results/caribbeanrisk/training/data.nc")
# ds = xr.open_dataset("/Users/alison/Documents/DPhil/github.nosync/hazGAN2/results/caribbeanrisk/processing/data_all.nc")
ds
ds['sro'] = ds['anomaly'].sel(field="sro")
ds['tp'] = ds['anomaly'].sel(field="tp")
ds['ws'] = ds['anomaly'].sel(field="ws")
# %%
fig, axs = plt.subplots(1, 3, figsize=(16, 4))

ds['tp'].isel(time=0).plot(cmap='Blues', ax=axs[0])
ds['sro'].isel(time=0).plot(cmap='Blues', ax=axs[1])
ds['ws'].isel(time=0).plot(cmap='Spectral_r', ax=axs[2])
# %%
from shapely.geometry import box
datadir = "/Users/alison/Documents/DPhil/data/hydrobasins/hybas_as_lev01-12_v1c"
import geopandas as gpd
# %%
LEVEL = 8
EXTENT = box(80, 10, 95, 25)
file = f"hybas_as_lev{str(LEVEL).zfill(2)}_v1c.shp"
gdf = gpd.read_file(os.path.join(datadir, file))
gdf = gdf.to_crs("EPSG:4326")
basins = gdf.clip(EXTENT)
# %%

import rioxarray

ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
ds.rio.write_crs("epsg:4326", inplace=True)

# %%
clipped = ds.rio.clip(gdf.geometry, gdf.crs, drop=True)
# %%
clipped['sro'].isel(time=0).plot(cmap='Blues')
# %%
runoff = clipped['sro'].to_dataframe().reset_index()[['time', 'lon', 'lat', 'sro']]
runoff = gpd.GeoDataFrame(runoff, geometry=gpd.points_from_xy(runoff.lon, runoff.lat))
runoff = runoff.set_crs("EPSG:4326")
# %%
runoff = runoff.dropna(subset='sro')
runoff.plot('sro')
# %%
gdf.plot("SORT")
# %%
intersected = gpd.overlay(basins, runoff, how="")
# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
runoff.plot(ax=ax, column='sro', cmap='Blues', legend=True)
basins.plot(ax=ax, color='none', edgecolor='black')
# intersected.plot(ax=ax, column='sro', cmap='Blues', legend=True)

# %%
# %%
intersection = gpd.sjoin(basins, runoff, how="left", predicate="contains")
# basins.sindex.valid_query_predicates
# %%
intersection = intersection.drop(columns=['index_right', 'lon', 'lat'])
intersection = intersection[["HYBAS_ID","NEXT_DOWN", "SUB_AREA", "SORT", "geometry", "sro"]].groupby(["HYBAS_ID","NEXT_DOWN", "SUB_AREA", "SORT", "geometry"]).sum()
intersection = intersection.reset_index()
# %%
intersection
# %%
intersection = gpd.GeoDataFrame(intersection, geometry=intersection.geometry)
# %%
intersection.plot("sro", cmap="Blues")
# %%


# %%
# 1 is most downstream
intersection = intersection.sort_values(by="SORT", ascending=False)

# %%
intersection["volume"] = intersection["sro"] * intersection["SUB_AREA"]
intersection.plot("volume", cmap="Blues", legend=True)


# %%
from matplotlib.colors import LogNorm
VOI = "sro" # "volume" or "sro"
fig, axs = plt.subplots(1, 2, figsize=(16, 8))

# propagate sro upstream --> downstream by summation
totals = intersection.copy().set_index("HYBAS_ID")
totals.sort_values(by=VOI, ascending=False).head()
totals.plot(VOI, cmap="Blues", legend=True, ax=axs[0])
flow_accumulator = totals.copy()

i = 0
while (len(flow_accumulator) > 0):
# while i < 3:
    print(i, len(flow_accumulator))
    flow_accumulator = flow_accumulator[["NEXT_DOWN", VOI]].groupby("NEXT_DOWN").agg(inflow=(VOI, "sum")).reset_index()
    flow_accumulator = flow_accumulator.rename(columns={"NEXT_DOWN": "HYBAS_ID"})
    flow_accumulator = flow_accumulator.set_index("HYBAS_ID")

    # add new inflows to basin volumes
    totals = totals.merge(flow_accumulator[["inflow"]], how="left", left_index=True, right_index=True)
    totals["inflow"] = totals["inflow"].fillna(0)
    totals[VOI] = totals[VOI] + totals["inflow"]
    totals = totals.drop(columns=["inflow"])

    # assign updated basin volumes to flow accumulatoe
    flow_accumulator = flow_accumulator.merge(totals[["NEXT_DOWN", VOI]], how="left", left_index=True, right_index=True)
    i += 1

totals.plot(VOI, cmap="Blues", norm=LogNorm(), legend=True, ax=axs[1])
# same plot but log colorscale
totals.sort_values(by=VOI, ascending=False).head()
# %%

# %%
