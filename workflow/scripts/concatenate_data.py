"""
Concatenate data from multiple netCDF files into a single netCDF file.
"""

import numpy as np
import xarray as xr
import time
import logging


logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

if __name__ == "__main__":
        start = time.time()
        logging.info("Starting resampling script.")

        # snakemake params
        INPUTS  = snakemake.input.netcdfs
        OUTPUT = snakemake.output.netcdf

        # load datasets
        ds = xr.open_mfdataset(INPUTS, chunks={"time": "500MB"}, engine="netcdf4")

        # add a grid variable for indexing
        h, w = ds.sizes['lat'], ds.sizes['lon']
        grid = np.arange(0, h * w, 1).reshape(h, w)
        grid = xr.DataArray(grid, dims=["lat", "lon"], coords={"lat": ds.lat[::-1], "lon": ds.lon})
        ds['grid'] = grid

        logging.info(f"Saving to {OUTPUT}")
        ds.to_netcdf(OUTPUT)

        logging.info(f"Saved {OUTPUT}")
        logging.info(f"Resampling took {time.time() - start} seconds.")