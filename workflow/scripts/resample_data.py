"""
Resample ERA5 (or other gridded) data to another lower resolution using GDAL
"""
# %%
import os
import numpy as np
import xarray as xr
from tqdm import tqdm
import argparse
import time
import logging
import glob


logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

if __name__ == "__main__":
    start = time.time()
    logging.info("Starting resampling script.")

    # snakemake params
    YEAR  = int(snakemake.params.year)
    RESX  = snakemake.params.resx
    RESY  = snakemake.params.resy
    INPUT  = snakemake.input.netcdf
    OUTPUT = snakemake.output.netcdf
    FIELDS = snakemake.params.fields

    # make output directory
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

    # find input files
    fields = list(FIELDS.keys())

    # load data for year
    logging.info(f"Resampling {YEAR} to {RESX}x{RESY}")
    ds_orig = xr.open_dataset(INPUT, engine="netcdf4")
    times = ds_orig.time.values

    # resample each field separately
    resampled_datasets = []
    for field, methods in FIELDS.items():
        # sample field
        method = methods["agg"]
        tmpfile = OUTPUT.replace('.nc', f'_{field}.nc')
        command = f'gdalwarp -t_srs EPSG:4326 -ts {RESX} {RESY} -r {method} -overwrite -of netCDF "NETCDF:{INPUT}:{field}" {tmpfile}'
        logging.info(f"Submitting GDAL command: {command}")
        os.system(command)

        # add resampled field to datasets holder
        ds_tmp = xr.open_dataset(tmpfile)
        bands = [var for var in ds_tmp.data_vars if 'Band' in var]
        ds_tmp = ds_tmp[bands].to_array('time', name=field).to_dataset().assign_coords(time=times)
        resampled_datasets.append(ds_tmp)

        # tidy
        os.remove(tmpfile)
            
    # merge datasets
    ds = xr.merge(resampled_datasets)
    ds.to_netcdf(OUTPUT)
    logging.info(f"Saved {OUTPUT}")
    logging.info(f"Resampling took {time.time() - start} seconds.")
#%%

