"""Make 64x64 plots of data

Note: this is before adding back monthly medians.
"""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from snakemake.script import snakemake

from src import plotting
from src import statistics

import ssl # fix cartopy issue
ssl._create_default_https_context = ssl._create_stdlib_context


if __name__ == "__main__":
    TRAIN = snakemake.input.train
    GENER = snakemake.input.generated
    FIELDS = snakemake.params.fields
    SHUFFLE = snakemake.params.shuffle
    OUTDIR  = snakemake.output.outdir

    os.makedirs(OUTDIR, exist_ok=True)

    # load data and samples
    train = xr.open_dataset(TRAIN)
    train_x = train["anomaly"].values
    train_u = train["uniform"].values
    train_g = statistics.gumbel(train_u)

    gener = xr.open_dataset(GENER)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values
    gener_g = statistics.gumbel(gener_u)

    if SHUFFLE:
        train_ids = np.random.permutation(train_u.shape[0])
        gener_ids = np.random.permutation(gener_u.shape[0])
    else:
        # train_ids = np.arange(train_u.shape[0])
        # gener_ids = np.arange(gener_u.shape[0])
        
        # sort by ascending wind speed
        # ! HARDCODED
        train_ids = np.argsort(np.max(train_u[..., 0], axis=(1, 2)))[::-1]
        gener_ids = np.argsort(np.max(gener_u[..., 0], axis=(1, 2)))[::-1]

    # make the samples
    for i, FIELD in enumerate(FIELDS.keys()):
        METRIC = FIELDS[FIELD]["units"]
        CMAP   = FIELDS[FIELD]["cmap"]

        figa = plotting.samples.plot(gener_g[gener_ids], train_g[train_ids], field=i, title="", cbar_label="", cmap=CMAP, ndecimals=0)
        figb = plotting.samples.plot(gener_u[gener_ids], train_u[train_ids], field=i, title="", cbar_label="", cmap=CMAP, ndecimals=1)
        figc = plotting.samples.plot(gener_x[gener_ids], train_x[train_ids], field=i, title="", cbar_label=METRIC, cmap=CMAP, alpha=1e-6);

        figa.savefig(os.path.join(OUTDIR, f"{FIELD}_gumbel.png"), dpi=300, bbox_inches="tight")
        figb.savefig(os.path.join(OUTDIR, f"{FIELD}_uniform.png"), dpi=300, bbox_inches="tight")
        figc.savefig(os.path.join(OUTDIR, f"{FIELD}_anomaly.png"), dpi=300, bbox_inches="tight")

        dates = {
            "generated": gener.time.values[gener_ids],
            "training": train.time.values[train_ids],
        }
        # save to .txt file
        with open(os.path.join(OUTDIR, f"{FIELD}_dates.txt"), "w") as f:
            for data in dates["training"]:
                f.write(f"Training date: {data}\n")


