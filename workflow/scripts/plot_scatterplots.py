# %% 
import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from shapely.geometry import Point

import logging

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


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
        )
    
    TRAIN = snakemake.input.train
    GENER = snakemake.input.generated
    DIR  = snakemake.output.outdir
    DO_SUBSET = snakemake.params.do_subset
    THRESH = snakemake.params.event_subset
    YMIN = snakemake.params.ymin
    YMAX = snakemake.params.ymax
    XMIN = snakemake.params.xmin
    XMAX = snakemake.params.xmax
    POIS = snakemake.params.pois
    CMAP = snakemake.params.cmap
    FIELD_LABELS = snakemake.params.channel_labels

    os.makedirs(DIR, exist_ok=True)

    # load data and samples
    train = xr.open_dataset(TRAIN)
    if DO_SUBSET:
        train['intensity'] = getattr(
            train.sel(field=THRESH["field"]).anomaly,
            THRESH["func"])(dim=['lon', 'lat']
                            )
        mask = (train['intensity'] > THRESH["value"]).values
        idx  = np.where(mask)[0]
        train   = train.isel(time=idx)
    
    train_x = train["anomaly"].values
    train_u = train["uniform"].values

    gener = xr.open_dataset(GENER)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values
    
    # match observation points to indices in the xarray data
    print(f"Observation points: {POIS}")
    ops = op2idx(POIS, train_x[0, ..., 0], extent=[XMIN, XMAX, YMIN, YMAX])
    pixels = [ops["chittagong"], ops["dhaka"]]

    FIELDS = FIELD_LABELS.keys()

    for i, FIELD in enumerate(FIELDS):
        fig = plotting.scatter.plot(gener_x, train_x, field=i, pixels=pixels, s=10,
                        cmap=CMAP, xlabel="Chittagong", ylabel="Dhaka")

        fig.suptitle(FIELD_LABELS[FIELD].capitalize(), y=1.05, fontsize=14, fontweight='bold')
        fig.savefig(os.path.join(DIR, f"field_{FIELD}.png"), dpi=300)
        plt.close(fig)

# %%
