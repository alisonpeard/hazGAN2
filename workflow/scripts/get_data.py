
"""
Load gridded data from hourly netcdf files, resample to daily aggregates, and save to
a single netcdf file in the target directory.
"""
import os
from pathlib import Path
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
    text = f"Data summary:\n\n"
    text += f"Dimensions: {data.dims}\n"
    text += f"Coordinates: {data.coords}\n"
    text += f"Size: {data.nbytes * 1e-6:.2f} MB\n"
    text += f"Variables: {data.data_vars}\n\n"
    logging.info(text)


def main(input, output, params):
    """Load data variables and process to daily aggregates."""
    logging.info("Starting data acquisition script.\n")
    start = time.time()
    dataset = getattr(datasets, params.dataset)

    input_files = dict()
    for field, field_meta in params.fields.items():
        for arg in field_meta["init"]["args"]:
            if arg not in input_files:
                input_file_pattern = dataset.get_input_file_pattern(input.indir, arg)
                arg_files = glob(input_file_pattern)
                arg_files = dataset.filter_files(
                    arg_files, params.year,
                    antecedent_buffer_days=params.antecedent_buffer_days
                )
                input_files[arg] = arg_files

    logging.debug(f"Input file {-1}: {input_files[arg][-1]}")

    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        counter = [0]
        def preprocess(ds, params=params):
            counter[0] += 1
            logging.info(f"Preprocessing file {counter[0]} of {len(input_files)}")
            # convert lon from 0-360 to -180 to 180
            if 'longitude' in ds.coords and ds.longitude.max() > 180:
                ds.coords['longitude'] = ((ds.coords['longitude']+180)%360-180)
                ds = ds.sortby(ds.longitude)
            ds = ds.sel(latitude=slice(params.ymax, params.ymin), longitude=slice(params.xmin, params.xmax))
            return ds
        
        ds_list = []
        for arg, arg_files in input_files.items():
            logging.info(f"Loading mfdataset for {arg}")
            ds = xr.open_mfdataset(
                arg_files,
                engine='cfgrib',
                preprocess=preprocess,
                combine="nested",
                concat_dim="valid_time",
                parallel=False,
                chunks="auto",
                backend_kwargs={
                    'time_dims': ('valid_time',),
                    'indexpath': ''
                }
            )
            ds_list.append(ds)
        
    data = xr.merge(ds_list).rename({'valid_time': 'time'})
    log_data_summary(data)

    logging.info("Computing data variables...")

    if params.antecedent_buffer_days:
        t0 = data["time"].min().data
        tn = data["time"].max().data

        logging.info(f"Data time range: {t0} to {tn}.")

        buffer = np.timedelta64(params.antecedent_buffer_days+1, 'D')
        data = data.sortby("time")
        t0 = np.datetime64(f"{params.year}-01-01")
        tn = np.datetime64(f"{params.year}-12-31")
        t0 -= buffer
        data = data.sel(time=slice(t0, tn))
        
        logging.info(f"Clipped data to time range: {t0} to {tn}.")
        
        t0 = data["time"].min().data
        tn = data["time"].max().data

        logging.info(f"Data time range: {t0} to {tn}.")

    logging.info(f"Loading parameters from {input.params}.\n")
    theta = xr.open_dataset(input.params)
    theta = preprocess(theta)

    logging.info("Resampling data to daily aggregates.\n")
    resampled = {}
    for field, config in params.fields.items():
        func       = config["init"]["func"]
        args       = config["init"]["args"]

        logging.info(f"Applying {field} = {func}{*args,}.")
        data[field] = getattr(funcs, func)(data, *args, params=theta)

    for field, config in params.fields.items():
        hfunc      = config["hfunc"]["func"]
        hfunc_args = config["hfunc"]["args"]

        logging.info(f"Resampling {field} = {hfunc}{*hfunc_args,}.")
        resampled[field] = data.resample(time='1D').apply(
            lambda grouped: getattr(funcs, hfunc)(grouped, *hfunc_args)
        )
    
    data_resampled = xr.Dataset(resampled)
    data_resampled = data_resampled.chunk({
        'time': 30,
        'latitude': -1, 
        'longitude': -1
    })

    if params.antecedent_buffer_days:
        t0 = np.datetime64(f"{params.year}-01-01")
        tn = np.datetime64(f"{params.year}-12-31")
        data_resampled = data_resampled.sel(time=slice(t0, tn))

    log_data_summary(data_resampled)

    # data export settings
    compression = {'zlib': True, 'complevel': 1}
    encoding = {var: compression for var in data_resampled.data_vars}

    # save data to netcdf
    logging.info(f"Saving data to {output.netcdf}\n")
    data_resampled.to_netcdf(output.netcdf, engine='netcdf4', encoding=encoding)
    logging.info(f"Saved. File size: {os.path.getsize(output.netcdf) * 1e-6:.2f} MB\n")

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
