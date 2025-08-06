import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from shapely.geometry import Point

import logging

from src import funcs
from src import plotting


def op2idx(ops:dict, data:np.ndarray, extent:list):
    """Convert observation points to indices in a 2D array."""
    h, w = data.shape

    lons = np.linspace(extent[0], extent[1], w)
    lats = np.linspace(extent[2], extent[3], h)

    coords = np.array([Point(lon, lat) for lon in lons for lat in lats])

    op_idx = {}
    for op, loc in ops.items():
        idx = np.argmin([coord.distance(Point(loc)) for coord in coords])
        op_idx[op] = idx

    return op_idx


def subset_func(ds:xr.Dataset, subset:dict):
    """Subset the dataset using function and threshold."""
    func = getattr(funcs, subset["func"])
    args = subset["args"]
    thresh = subset["value"]

    logging.info(f"Subsetting {func}{*args,} with threshold {thresh}.")

    for arg in args:
        ds[arg] = ds.sel(field=arg).anomaly

    intensity = func(ds, *args, dim=["lon", "lat"])

    for arg in args:
        ds = ds.drop_vars(arg)

    mask = (intensity > thresh).values
    idx = np.where(mask)[0]
    return ds.isel(time=idx)


def main(input, output, params):

    os.makedirs(output.outdir, exist_ok=True)

    # load data and samples
    train = xr.open_dataset(input.train)
    if params["subset"]["do"]:
        # subset train by threshold
        train = subset_func(train, params["subset"])
        logging.info(f"\nExtracted {train.time.size} images from train.")
    
    train_x = train["anomaly"].values
    train_u = train["uniform"].values

    gener = xr.open_dataset(input.generated)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values
    
    # match observation points to indices in the xarray data
    print(f"Observation points: {params.pois}")
    ops = op2idx(
        params.pois, train_x[0, ..., 0],
        extent=[params.xmin, params.xmax, params.ymin, params.ymax]
        )
    pixels = [ops["chittagong"], ops["dhaka"]]

    _fields = params.channel_labels.keys()

    for i, _field in enumerate(_fields):
        fig = plotting.scatter.plot(gener_x, train_x, field=i, pixels=pixels, s=10,
                        cmap=params.cmap, xlabel="Chittagong", ylabel="Dhaka")

        fig.suptitle(params.channel_labels[_field].capitalize(), y=1.05, fontsize=14, fontweight='bold')
        fig.savefig(os.path.join(output.outdir, f"field_{_field}.png"), dpi=300)
        plt.close(fig)


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
        )
    
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params

    main(input, output, params)

    logging.info("Plotting completed successfully.")