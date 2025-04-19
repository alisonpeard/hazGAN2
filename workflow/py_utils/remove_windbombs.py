"""
This script removes ERA5 "wind bombs", see ERA5 documentation for more details.
This is a little hardcoded, so might need to remove this for other projects. Including for completeness.
"""

import numpy as np
import xarray as xr
import time
import logging


def rescale(x:np.ndarray) -> np.ndarray:
    return (x - x.min() / (x.max() - x.min()))

def rescale_vector(x:np.ndarray) -> np.ndarray:
    return (x - x.min(axis=(1, 2), keepdims=True)) / (x.max(axis=(1, 2), keepdims=True) - x.min(axis=(1, 2), keepdims=True))

def frobenius_vector(test:np.ndarray, template:np.ndarray) -> np.ndarray:
    sum_ = np.sum(template * test, axis=(1, 2))
    norms = np.linalg.norm(template) * np.linalg.norm(test, axis=(1, 2))
    similarities = sum_ / norms
    return similarities

def get_similarities(ds:xr.Dataset, template:np.ndarray) -> np.ndarray:
    """Get similarities between a template and dataset."""
    template = rescale(template)
    tensor = ds['ws'].data
    tensor = rescale_vector(tensor)
    similarities = frobenius_vector(tensor, template)
    
    return similarities # np.array(similarities)


logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

if __name__ == "__main__":
        start = time.time()
        logging.info("Starting resampling script.")

        # snakemake parameters
        INPUT = snakemake.input.netcdf
        OUTPUT = snakemake.output.netcdf
        WINDBOMB = snakemake.output.windbomb
        THRESHOLD = snakemake.params.threshold

        # load the data
        ds = xr.open_dataset(INPUT, chunks={"time": "500MB"}, engine="netcdf4")
        logging.info("Loaded data.")

        # find the first wind bomb--usually corresponds to highest wind speed
        ds['maxwind'] = ds['ws'].max(dim=['lat', 'lon'])
        ds = ds.sortby("maxwind", ascending=False)
        template = ds.isel(time=0).ws.data
        logging.info("Calculated first wind bomb.")

        # calculate similatities
        similarities = get_similarities(ds, template).compute()
        nbombs = sum(similarities >= THRESHOLD)
        mask = (similarities <= THRESHOLD)
        ds_filtered = ds.isel(time=mask)
        ds_filtered = ds_filtered.sortby('time', ascending=True)

        # save outputs
        np.save(WINDBOMB, template.compute())
        ds_filtered.to_netcdf(OUTPUT, mode="w")