"""Make 64x64 plots of data

Note: this is before adding back monthly medians.
"""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from snakemake.script import snakemake

from hazGAN.plotting import samples
from hazGAN.statistics import gumbel

import ssl # fix cartopy issue
ssl._create_default_https_context = ssl._create_stdlib_context


if __name__ == "__main__":
    TRAIN = snakemake.input.train
    FIELDS = snakemake.params.fields
    SHUFFLE = snakemake.params.shuffle
    OUTDIR  = snakemake.output.outdir

    os.makedirs(OUTDIR, exist_ok=True)

    # load data and samples
    train = xr.open_dataset(TRAIN)
    train_x = train["anomaly"].values
    train_u = train["uniform"].values
    train_g = gumbel(train_u)

    if SHUFFLE:
        train_ids = np.random.permutation(train_u.shape[0])
    else:
        train_ids = np.argsort(np.max(train_u[..., 0], axis=(1, 2)))[::-1]

    # make the samples
    for i, FIELD in enumerate(FIELDS.keys()):
        METRIC = FIELDS[FIELD]["units"]
        CMAP   = FIELDS[FIELD]["cmap"]

        figa = samples.plot(train_g[train_ids], train_g[train_ids], field=i, title="", cbar_label="", cmap=CMAP, ndecimals=0)
        figb = samples.plot(train_u[train_ids], train_u[train_ids], field=i, title="", cbar_label="", cmap=CMAP, ndecimals=1)
        figc = samples.plot(train_x[train_ids], train_x[train_ids], field=i, title="", cbar_label=METRIC, cmap=CMAP, alpha=1e-6);

        figa.savefig(os.path.join(OUTDIR, f"{FIELD}_gumbel.png"), dpi=300, bbox_inches="tight")
        figb.savefig(os.path.join(OUTDIR, f"{FIELD}_uniform.png"), dpi=300, bbox_inches="tight")
        figc.savefig(os.path.join(OUTDIR, f"{FIELD}_anomaly.png"), dpi=300, bbox_inches="tight")

        dates = {
            "training": train.time.values[train_ids],
        }
        # save to .txt file
        with open(os.path.join(OUTDIR, f"{FIELD}_dates.txt"), "w") as f:
            for data in dates["training"]:
                f.write(f"Training date: {data}\n")


