
"""
Load ERA5 data from hourly netcdf files, resample to daily aggregates, and save to a single netcdf file in the target directory.
"""
#%%
import os
import sys
from environs import Env
import dask
import time
from glob import glob
from pprint import pp as prettyprint
import xarray as xr
import numpy as np
import logging


logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

if __name__ == '__main__':
    start = time.time()
    logging.info("Starting data acquisition script.")

    # check we are on the cluster
    env = snakemake.config.get("env", os.environ.get("env", "local"))
    if env != "cluster":
        logging.error(f"Error: This script can only be run in the cluster environment. Current env: {env}", file=sys.stderr)
        sys.exit(1)

    # snakemake params
    YEAR       = int(snakemake.params.year)
    INDIR      = snakemake.params.indir
    XMIN       = snakemake.params.xmin
    XMAX       = snakemake.params.xmax
    YMIN       = snakemake.params.ymin
    YMAX       = snakemake.params.ymax
    FIELDS     = snakemake.params.fields
    FIELD_LONG = snakemake.params.field_long
    OUTPUT     = snakemake.output.netcdf

    # load data for year
    files = []
    for field in FIELD_LONG.values():
        field_files = glob(os.path.join(INDIR, field, 'nc', '*'))
        field_files = [f for f in field_files if str(YEAR) in f]
        files += field_files
    logging.info(f"Found {len(files)} files for {YEAR}...")

    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        data = xr.open_mfdataset(files, engine='netcdf4', chunks={"time": "500MB", 'longitude': '500MB', 'latitude': '500MB'})
        data = data.sel(longitude=slice(XMIN, XMAX), latitude=slice(YMAX, YMIN))
    logging.info("Data loaded.")

    # resample data to daily
    resampled = {}
    for var, func in FIELDS.items():
        resampled[var] = getattr(data[var].resample(time='1D'), func)()
        logging.info(f"Resampled {var} using {func}.")
    data_resampled = xr.Dataset(resampled)

    # save data
    chunk_size = {'time': '50MB'}
    data_resampled = data_resampled.chunk(chunk_size)
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    logging.info(f"Saving data to {OUTPUT}")
    data_resampled.to_netcdf(OUTPUT, engine='netcdf4')
    logging.info(f"Saved. File size: {os.path.getsize(OUTPUT) * 1e-6:.2f} MB")

    # tidy up
    data.close()
    data_resampled.close()

    # Verify time encoding
    data = xr.open_dataset(OUTPUT, decode_times=False, engine='netcdf4')
    times = data.isel(time=slice(0, 4)).time.data
    logging.info(f"Time encoding: {times[0]}, {times[1]}, ...")
    logging.info(f"Encoding metadata: {data.time.attrs}")
    data.close()

    end = time.time()
    logging.info(f"Time taken: {end - start:.2f} seconds.")
# %%
