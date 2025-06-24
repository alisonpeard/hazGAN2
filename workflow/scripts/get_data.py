
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



def main(input, output, params):
    start = time.time()
    logging.info("Starting data acquisition script.")
    # load dataset module
    dataset = import_module(f"py_utils.datasets.{params.dataset}")

    # load and clip to area of interest
    files = []
    for field, info in params.fields.items():
        infields = info.get("args", [])
        for infield in infields:
            pattern = dataset.parse_input_pattern(input.indir, infield)
            field_files = glob(pattern)
            field_files = dataset.filter_files(field_files, params.year)
            files += field_files
    
    logging.info(f"Found {len(files)} files for {params.year}.")
    for i, file in enumerate(files):
        logging.info(f"File {i}: {file}")

    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
        # data = xr.open_mfdataset(files, engine='netcdf4',
        #                          chunks={
        #                              "time": "500MB",
        #                              params.timecol: "500MB",
        #                              'longitude': '500MB',
        #                              'latitude': '500MB'
        #                              })
        datasets = []
        for file in files:
            data_i = xr.open_dataset(file, engine='netcdf4', chunks={
                "time": "500MB",
                params.timecol: "500MB",
                'longitude': '500MB',
                'latitude': '500MB'
            })

            if params.timecol in data_i.coords:
                logging.warning(f"Removing coordinate {params.timecol} from dataset {file}.")
                data_i = data_i.rename({params.timecol: "time"})

            data_i = dataset.clip_to_bbox(data_i, params.xmin, params.xmax, params.ymin, params.ymax)
            
            logging.info(f"\nData loaded with {len(data_i.time)} time steps, {len(data_i.latitude)} latitudes, and {len(data_i.longitude)} longitudes.")
            logging.info(f"Dimensions: {data_i.dims} and variables: {data_i.data_vars}")

            datasets.append(data_i)

        data = xr.merge(datasets)
        logging.info(f"Data loaded with {len(data.time)} time steps, {len(data.latitude)} latitudes, and {len(data.longitude)} longitudes.")
        logging.info(f"Dimensions: {data.dims} and variables: {data.data_vars}")
        
    logging.info("\n\nData loaded.")

    # log data summary
    logging.info("Data summary:")
    logging.info(f"Data variables: {data.data_vars}")
    logging.info(f"Data coordinates: {data.coords}")
    logging.info(f"Data dimensions: {data.dims}")

    logging.info(f"Loading parameters from {input.params}.")
    theta = xr.open_dataset(input.params)
    theta = theta.sel(longitude=slice(params.xmin, params.xmax), latitude=slice(params.ymax, params.ymin))

    logging.info("Resampling data to daily aggregates.")
    resampled = {}
    for field, config in params.fields.items():
        logging.info(f"Processing {field}.")
        infields = config.get("args", [])
        func     = config.get("func", "identity")
        hfunc    = config.get("hfunc", "mean")

        logging.info(f"Applying {field} = {func}{*infields,}.")
        data[field] = getattr(funcs, func)(*[data[i] for i in infields], params=theta)

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
    os.makedirs(os.path.dirname(output.netcdf), exist_ok=True)
    logging.info(f"Saving data to {output.netcdf}")
    data_resampled.to_netcdf(output.netcdf, engine='netcdf4', encoding=encoding)
    logging.info(f"Saved. File size: {os.path.getsize(output.netcdf) * 1e-6:.2f} MB")

    # print data summary to logs
    logging.info("Data summary:")
    logging.info(f"Data variables: {data_resampled.data_vars}")
    logging.info(f"Data coordinates: {data_resampled.coords}")

    # tidy up
    data.close()
    data_resampled.close()

    # verify time encoding
    data = xr.open_dataset(output.netcdf, decode_times=False, engine='netcdf4')
    times = data.isel(time=slice(0, 4)).time.data
    logging.info(f"Time encoding: {times[0]}, {times[1]}, ...")
    logging.info(f"Encoding metadata: {data.time.attrs}")
    data.close()

    end = time.time()
    logging.info(f"Time taken: {end - start:.2f} seconds.")


if __name__ == '__main__':
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    main(input, output, params)
