import os
import yaml
import geopandas as gpd
import xarray as xr
from shapely.geometry import box
from hazGAN.mangrove_demo import mangroveDamageModel

MANGROVES = "/Users/alison/Documents/DPhil/data/gmw-v3-2020.nosync/gmw_v3_2020_vec.gpkg"

if __name__ == "__main__":
    wd = os.path.join("..", "results", "mangroves")
    os.makedirs(wd, exist_ok=True)

    with open(os.path.join("..", "config.yaml"), 'r') as stream:
        config = yaml.safe_load(stream)

    if os.path.exists(os.path.join(wd, "mangrove_grid.nc")):
        print("mangrove_grid.nc already exists, skipping mangrove grid generation.")
        exit(0)
    else:
        bay_of_bengal_crs = config['local_crs']
        xmin = config["longitude"]["min"]
        xmax = config["longitude"]["max"]
        ymin = config["latitude"]["min"]
        ymax = config["latitude"]["max"]

        # load reference data
        aoi = box(xmin, ymin, xmax, ymax)
        train_damages = xr.open_dataset(os.path.join(wd, "train_damages.nc"))
        grid_damages = train_damages.isel(time=0)

        #  clip mangroves to bay of bengal
        mangroves = gpd.read_file(MANGROVES, mask=aoi)
        mangroves = mangroves.set_crs(epsg=4326).drop(columns='PXLVAL')
        mangroves['area']  = mangroves.to_crs(bay_of_bengal_crs).area
        mangrove_centroids = mangroves.set_geometry(mangroves.centroid)
        mangrove_centroids = mangrove_centroids.sort_values(by='area', ascending=False)

        # intersect mangroves with grid
        mangroves = gpd.read_file(os.path.join("..", "resources", "mangroves.geojson"))
        model = mangroveDamageModel()
        mangrove_grid = model.intersect_mangroves_with_grid(mangroves, grid_damages)
        mangrove_grid.to_netcdf(os.path.join(wd, "mangrove_grid.nc"))