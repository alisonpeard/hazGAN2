
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

            ds = dataset.preprocess(ds)

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
    
    logging.info("Computing data variables...")
    log_data_summary(data)

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
    data_resampled.compute().to_netcdf(output.netcdf, engine='netcdf4', encoding=encoding)
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
