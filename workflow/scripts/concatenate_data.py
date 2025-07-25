"""
Concatenate data from multiple netCDF files into a single netCDF file.
"""

import numpy as np
import xarray as xr
import time
import logging


logging.basicConfig(
    filename=snakemake.log.file,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def main(input, output, params):
        start = time.time()
        logging.info("Starting resampling script.")

        ds = xr.open_mfdataset(input.netcdfs, chunks={"time": "500MB"}, engine="netcdf4")

        logging.info(f"Saving to {output.netcdf}")
        encoding = {
            var: {"zlib": True, "complevel": 5} for var in ds.data_vars
        }

        # remove any dates that are in the exclude list
        if params.exclude:
            exclude_dates = [np.datetime64(date) for date in params.exclude]
            date_mask = np.isin(ds.time.values, exclude_dates, invert=True)
            ds = ds.isel(time=date_mask)

        ds.to_netcdf(output.netcdf, encoding=encoding)

        logging.info(f"Saved {output.netcdf}")
        logging.info(f"Resampling took {time.time() - start} seconds.")


if __name__ == "__main__":
        input = snakemake.input
        output = snakemake.output
        params = snakemake.params
        main(input, output, params)
