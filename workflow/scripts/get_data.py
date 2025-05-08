
"""
Load ERA5 data from hourly netcdf files, resample to daily aggregates, and save to a single netcdf file in the target directory.
"""
import os
from glob import glob
import dask
import time
import xarray as xr
import logging

from importlib import import_module
from py_utils import funcs

logging.basicConfig(
    filename=snakemake.log.file, # try withou the [0] in case it's makin it a dir
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logging.info("Starting data acquisition script (before main).")

if __name__ == '__main__':
    start = time.time()
    logging.info("Starting data acquisition script.")

    # snakemake params
    INDIR      = snakemake.input.indir
    INPUT      = snakemake.input.params

    YEAR       = int(snakemake.params.year)
    XMIN       = snakemake.params.xmin
    XMAX       = snakemake.params.xmax
    YMIN       = snakemake.params.ymin
    YMAX       = snakemake.params.ymax
    FIELDS     = snakemake.params.fields
    TIMECOL    = snakemake.params.timecol
    DATASET    = snakemake.params.dataset

    OUTPUT     = snakemake.output.netcdf

    # load dataset module
    dataset = import_module(f"py_utils.datasets.{DATASET}")

    # load params and clip to XMIN, XMAX, YMIN, YMAX
    logging.info(f"Loading parameters from {INPUT}.")
    params = xr.open_dataset(INPUT)
    params = params.sel(longitude=slice(XMIN, XMAX), latitude=slice(YMAX, YMIN))

    # load data for year
    files = []
    for field, info in FIELDS.items():
        infields = info.get("args", [])
        for infield in infields:
            # find the input file(s)
            pattern = dataset.parse_input_pattern(INDIR, infield)
            field_files = glob(pattern)
            field_files = dataset.filter_files(field_files, YEAR)
            files += field_files
    
    logging.info(f"Found {len(files)} files for {YEAR}.")
    for i, file in enumerate(files):
        logging.info(f"File {i}: {file}")

    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        data = xr.open_mfdataset(files, engine='netcdf4',
                                 chunks={
                                     TIMECOL: "500MB",
                                     'longitude': '500MB',
                                     'latitude': '500MB'
                                     })
        data = dataset.clip_to_bbox(data, XMIN, XMAX, YMIN, YMAX)
        data = data.rename({TIMECOL: "time"})
    logging.info("Data loaded.")

    # log data summary
    logging.info("Data summary:")
    logging.info(f"Data variables: {data.data_vars}")
    logging.info(f"Data coordinates: {data.coords}")
    logging.info(f"Data dimensions: {data.dims}")

    # resample data to daily
    logging.info("Resampling data to daily aggregates.")
    resampled = {}
    for field, config in FIELDS.items():
        logging.info(f"Processing {field}.")
        infields = config.get("args", [])
        func     = config.get("func", "identity")
        hfunc    = config.get("hfunc", "mean")

        logging.info(f"Applying {field} = {func}{*infields,}.")
        data[field] = getattr(funcs, func)(*[data[i] for i in infields], params=params)

        logging.info(f"Finished, resampling as {hfunc}({field}).")
        resampled[field] = getattr(data[field].resample(time='1D'), hfunc)()
        
        logging.info(f"Finished {field}.")
    
    data_resampled = xr.Dataset(resampled)

    # log data summary
    logging.info("Resampled data summary:")
    logging.info(f"Data variables: {data_resampled.data_vars}")
    logging.info(f"Data coordinates: {data_resampled.coords}")
    logging.info(f"Data dimensions: {data_resampled.dims}")
    logging.info(f"Data size: {data_resampled.nbytes * 1e-6:.2f} MB")

    # data export settings
    chunk_size  = {'time': '50MB'}
    compression = {'zlib': True, 'complevel': 5}
    encoding = {var: compression for var in data_resampled.data_vars}

    # save data to netcdf
    data_resampled = data_resampled.chunk(chunk_size)
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    logging.info(f"Saving data to {OUTPUT}")
    data_resampled.to_netcdf(OUTPUT, engine='netcdf4', encoding=encoding)
    logging.info(f"Saved. File size: {os.path.getsize(OUTPUT) * 1e-6:.2f} MB")

    # print data summary to logs
    logging.info("Data summary:")
    logging.info(f"Data variables: {data_resampled.data_vars}")
    logging.info(f"Data coordinates: {data_resampled.coords}")

    # tidy up
    data.close()
    data_resampled.close()

    # verify time encoding
    data = xr.open_dataset(OUTPUT, decode_times=False, engine='netcdf4')
    times = data.isel(time=slice(0, 4)).time.data
    logging.info(f"Time encoding: {times[0]}, {times[1]}, ...")
    logging.info(f"Encoding metadata: {data.time.attrs}")
    data.close()

    end = time.time()
    logging.info(f"Time taken: {end - start:.2f} seconds.")
# %%
