
"""
Load gridded data from hourly netcdf files, resample to daily aggregates, and save to
a single netcdf file in the target directory.
"""
import os
import numpy as np
import sys
from glob import glob
import dask
import time
import xarray as xr
import logging
from src import funcs
from src import datasets


VERIFY_TIME_ENCODING = True  # enable for dev only


def log_data_summary(data):
    """Log a summary of the data."""
    logging.info("\n\nData summary:")
    logging.info(f"Data dimensions: {data.dims}")
    logging.info(f"Data coordinates: {data.coords}")
    logging.info(f"Data size: {data.nbytes * 1e-6:.2f} MB")
    logging.info(f"Data variables: {data.data_vars}\n")


def main(input, output, params):
    """Load data variables and process to daily aggregates."""
    logging.info("Starting data acquisition script.")
    start = time.time()
    dataset = getattr(datasets, params.dataset)

    input_files = set()
    for field, field_meta in params.fields.items():
        args = field_meta["init"]["args"]
        for arg in args:
            input_file_pattern = dataset.get_input_file_pattern(input.indir, arg)
            arg_files = glob(input_file_pattern)
            arg_files = dataset.filter_files(arg_files, params.year, antecedent_buffer_days=params.antecedent_buffer_days)
            input_files.update(arg_files)

    for i, file in enumerate(input_files):
        logging.info(f"Input file {i}: {file}")

    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        def preprocess(ds, params=params):
            """Rename time coordinate if necessary."""
            if params.xmin < 0:
                ds = funcs.convert_360_to_180(ds)
            ds = dataset.clip_to_bbox(
                ds, params.xmin, params.xmax, params.ymin, params.ymax
            )
            if params.timecol in ds.coords and params.timecol != "time":
                ds = ds.rename({params.timecol: "time"})
            if params.antecedent_buffer_days:
                if "time" in ds.coords:
                    buffer = np.timedelta64(params.antecedent_buffer_days, 'D')
                    ds = ds.sortby("time")
                    t0 = np.datetime64(f"{params.year}-01-01")
                    tn = np.datetime64(f"{params.year}-12-31")
                    t0 -= buffer
                    ds = ds.sel(time=slice(t0, tn))
                    logging.info(f"Clipped data to time range: {t0} to {tn}.")
            return ds

        data = xr.open_mfdataset(
            input_files,
            engine='netcdf4',
            preprocess=preprocess,
            chunks={
                "time": "500MB",
                'longitude': '500MB',
                'latitude': '500MB'
                })
    log_data_summary(data)

    # select expver=1 is expver in coords
    if "expver" in data.coords:
        logging.info("\nSelecting expver=1 (ERA5).")
        data = data.drop_vars('expver')
        log_data_summary(data)

    logging.info(f"\nLoading parameters from {input.params}.")
    theta = xr.open_dataset(input.params)
    theta = preprocess(theta)

    logging.info("\nResampling data to daily aggregates.")
    resampled = {}
    for field, config in params.fields.items():
        func       = config["init"]["func"]
        args       = config["init"]["args"]
        hfunc      = config["hfunc"]["func"]
        hfunc_args = config["hfunc"]["args"]

        logging.info(f"Applying {field} = {func}{*args,}.")
        data[field] = getattr(funcs, func)(data, *args, params=theta)

        logging.info(f"Resampling {field} = {hfunc}{*hfunc_args,}.")
        resampled[field] = data.resample(time='1D').apply(
            lambda grouped: getattr(funcs, hfunc)(grouped, *hfunc_args)
        )
    
    data_resampled = xr.Dataset(resampled)
    log_data_summary(data_resampled)

    # data export settings
    compression = {'zlib': True, 'complevel': 1}
    encoding = {var: compression for var in data_resampled.data_vars}

    # save data to netcdf
    logging.info(f"\n\nSaving data to {output.netcdf}")
    data_resampled.to_netcdf(output.netcdf, engine='netcdf4', encoding=encoding)
    logging.info(f"Saved. File size: {os.path.getsize(output.netcdf) * 1e-6:.2f} MB")

    data.close()
    data_resampled.close()

    if VERIFY_TIME_ENCODING:
        data = xr.open_dataset(output.netcdf, decode_times=False, engine='netcdf4')
        times = data.isel(time=slice(0, 4)).time.data
        logging.info(f"Time encoding: {times[0]}, {times[1]}, ...")
        logging.info(f"Encoding metadata: {data.time.attrs}")
        data.close()

    end = time.time()
    logging.info(f"Time taken: {end - start:.2f} seconds.")


if __name__ == '__main__':
    logging.basicConfig(
        # filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
        logging.FileHandler(snakemake.log.file),
        logging.StreamHandler(sys.stdout)
    ]
    )

    # process snakemake
    input  = snakemake.input
    output = snakemake.output
    params = snakemake.params

    main(input, output, params)
