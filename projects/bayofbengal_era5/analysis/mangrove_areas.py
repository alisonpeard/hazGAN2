"""Import mangrove data and clip to Bay of Bengal AOI.

Get area of each mangrove patch according to local CRS.
For each mangrove patch, the relative area of each covering grid 
cell is calculated 

I could probably add a small snakefile here too.
"""
# %%
import os
import yaml
import geopandas as gpd
import xarray as xr
import numpy as np
from shapely.geometry import box

import xagg as xa # Kalianpur 1975
import matplotlib.pyplot as plt

MANGROVES = "/Users/alison/Documents/DPhil/data/gmw-v3-2020.nosync/gmw_v3_2020_vec.gpkg"

if __name__ == "__main__":
    wd = os.path.join("..", "results")
    os.makedirs(wd, exist_ok=True)

    with open(os.path.join("..", "config.yaml"), 'r') as stream:
        config = yaml.safe_load(stream)

    if os.path.exists(os.path.join(wd, "mangroves", "mangrove_grid.nc")):
        print("mangrove_grid.nc already exists, skipping mangrove grid generation.")
        exit(0)
    else:
        crs  = config['local_crs']
        xmin = config["longitude"]["min"]
        xmax = config["longitude"]["max"]
        ymin = config["latitude"]["min"]
        ymax = config["latitude"]["max"]

        # load reference data
        aoi = box(xmin, ymin, xmax, ymax)
        grid_damages = xr.open_dataset(os.path.join(wd, "generated", "netcdf", "train.nc"))
        grid = grid_damages.isel(time=0, field=0)
        grid = grid.rio.write_crs("4326")
        grid = grid.rio.set_spatial_dims(x_dim='lon', y_dim='lat')

        # load and clip mangroves to RoI
        mangroves = gpd.read_file(MANGROVES, mask=aoi)
        mangroves = mangroves.set_crs(epsg=4326).drop(columns='PXLVAL')
        mangroves['area']  = mangroves.to_crs(crs).area
        
        # initialise plot
        fig, axs = plt.subplots(1, 2, figsize=(10, 10))
        ax = axs[0]
        grid.anomaly.plot(ax=ax, cmap="YlOrRd", add_colorbar=False)
        mangroves.plot("area", ax=ax, cmap="Greens", alpha=0.5, edgecolor='black', linewidth=0.1)
        # fig.savefig(os.path.join(wd, "figures", "mangroves", "input_data.pdf"), bbox_inches='tight')
        
        # intersect mangroves with grid
        weightmap = xa.pixel_overlaps(grid, mangroves)
        mangroves_gridded = weightmap.agg
        mangroves_gridded['npix'] = mangroves_gridded['pix_idxs'].apply(len)    
        mangroves_gridded['rel_area'] = mangroves_gridded['rel_area'].apply(lambda x: np.squeeze(x, axis=0))
        mangroves_gridded = mangroves_gridded.explode(['rel_area', 'pix_idxs'])

        # sum all relative mangrove areas in the same pixel
        mangroves_gridded['area'] = mangroves_gridded['area'] * mangroves_gridded['rel_area']
        mangroves_gridded = mangroves_gridded.groupby('pix_idxs').agg({'area': 'sum', 'coords': 'first'})

        # convert pd.DataFrame to xarray.Dataset
        lons = weightmap.source_grid['lon'].values
        lats = weightmap.source_grid['lat'].values
        mangroves_gridded = mangroves_gridded.reset_index()
        mangroves_gridded['lon'] = mangroves_gridded['pix_idxs'].apply(lambda j: lons[j])
        mangroves_gridded['lat'] = mangroves_gridded['pix_idxs'].apply(lambda i: lats[i])
        mangroves_gridded['lon'] = mangroves_gridded['lon'].astype(float)
        mangroves_gridded['lat'] = mangroves_gridded['lat'].astype(float)
        mangroves_gridded['area'] = mangroves_gridded['area'].astype(float)
        mangroves_gridded['area'] = mangroves_gridded['area'] * 1e-6 # convert to sqkm
        mangroves_gridded = mangroves_gridded.set_index(['lat', 'lon'])[['area']]
        mangroves_gridded = xr.Dataset.from_dataframe(mangroves_gridded)

        # process result
        mangroves_gridded = mangroves_gridded.rio.write_crs(crs)
        mangrove_grid = mangroves_gridded.rio.set_spatial_dims(x_dim='lon', y_dim='lat')

        # add to grid dataset
        zeros = 0 * grid.anomaly
        grid["mangrove_area"] = mangroves_gridded['area'] + zeros
        grid["mangrove_area"] = grid["mangrove_area"].fillna(0)

        ax = axs[1]
        grid.mangrove_area.plot(cmap="Greens", ax=ax, add_colorbar=False)
        mangroves.boundary.plot(ax=ax, color="k", linewidth=0.1)
        axs[0].set_title("Mangrove area in Bay of Bengal")
        axs[1].set_title("Mangrove area in Bay of Bengal (gridded)")
        fig.savefig(os.path.join(wd, "figures", "mangroves", "mangrove_area.pdf"), bbox_inches='tight')

        grid = grid.drop_vars(["anomaly", "gumbel", "uniform"])
        grid.to_netcdf(os.path.join(wd, "mangroves", "mangrove_areas.nc"))
        # %%