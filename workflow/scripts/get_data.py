
"""
Load ERA5 data from hourly netcdf files, resample to daily aggregates, and save to a single netcdf file in the target directory.
"""
import os
from glob import glob
import dask
import time
import xarray as xr
import logging
import py_utils.era5_utils as era5

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

    OUTPUT     = snakemake.output.netcdf

    # load params and clip to XMIN, XMAX, YMIN, YMAX
    logging.info(f"Loading parameters from {INPUT}.")
    params = xr.open_dataset(INPUT)
    params = params.sel(longitude=slice(XMIN, XMAX), latitude=slice(YMAX, YMIN))

    # load data for year
    files = []
    for field, info in FIELDS.items():
        infields = info.get("args", [])
        for infield in infields:
            field_files = glob(os.path.join(INDIR, era5.long_names[infield], 'nc', '*'))
            field_files = [f for f in field_files if str(YEAR) in f]
            files += field_files
    logging.info(f"Found {len(files)} files for {YEAR}.")

    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        data = xr.open_mfdataset(files, engine='netcdf4',
                                 chunks={
                                     TIMECOL: "500MB",
                                     'longitude': '500MB',
                                     'latitude': '500MB'
                                     })
        data = data.sel(longitude=slice(XMIN, XMAX), latitude=slice(YMAX, YMIN))
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
    for field, info in FIELDS.items():
        logging.info(f"Processing {field}.")
        infields = info.get("args", [])
        func = info.get("func", "identity")
        agg = info.get("agg", "mean")
        #! This is where output variables are created
        #! Add parameter loading here
        #! Untested
        logging.info(f"Applying {field} = {func}{*infields,}.")
        data[field] = getattr(era5, func)(*[data[i] for i in infields], params=params)
        logging.info(f"Finished, resampling {field}.")
        resampled[field] = getattr(data[field].resample(time='1D'), agg)()
        logging.info(f"Resampled using {agg}){field}).")
    data_resampled = xr.Dataset(resampled)

    # save data
    chunk_size = {'time': '50MB'}
    data_resampled = data_resampled.chunk(chunk_size)
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    logging.info(f"Saving data to {OUTPUT}")
    data_resampled.to_netcdf(OUTPUT, engine='netcdf4')
    logging.info(f"Saved. File size: {os.path.getsize(OUTPUT) * 1e-6:.2f} MB")

    # print data summary to logs
    logging.info("Data summary:")
    logging.info(f"Data variables: {data_resampled.data_vars}")
    logging.info(f"Data coordinates: {data_resampled.coords}")

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
