"""Make 64x64 plots of data"""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from snakemake.scripts import snakemake
from hazGAN.plotting import samples
from hazGAN.statistics import gumbel


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'


if __name__ == "__main__":
    TRAIN = snakemake.input.train
    GENER = snakemake.input.generated
    FIELDS = snakemake.params.fields
    SHUFFLE = snakemake.params.shuffle
    OUTDIR  = snakemake.output.outdir

    # load data and samples
    train = xr.open_dataset(TRAIN)
    train_x = train["anomaly"].values
    train_u = train["uniform"].values
    train_g = gumbel(train_x)

    gener = xr.open_dataset(GENER)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values
    gener_g = gumbel(gener_x)

    if SHUFFLE:
        train_ids = np.random.permutation(train_u.shape[0])
        gener_ids = np.random.permutation(gener_u.shape[0])
    else:
        train_ids = np.arange(train_u.shape[0])
        gener_ids = np.arange(gener_u.shape[0])

    # make the samples
    for i, FIELD in enumerate(FIELDS.keys()):
        METRIC = FIELDS[FIELD]["units"]
        CMAP   = FIELDS[FIELD]["cmap"]

        figa = samples.plot(gener_g[gener_ids], train_g[train_ids], field=i, title="", cmap=CMAP, ndecimals=0)
        figb = samples.plot(gener_u[gener_ids], train_u[train_ids], field=i, title="", cbar_label="", cmap=CMAP, ndecimals=1)
        figc = samples.plot(gener_x[gener_ids], train_x[train_ids], field=i, title="", cbar_label=METRIC, cmap=CMAP, alpha=1e-6);

        figa.savefig(os.path.join(OUTDIR, f"{FIELD}_gumbel.png", dpi=300))
        figb.savefig(os.path.join(OUTDIR, f"{FIELD}_uniform.png", dpi=300))
        figc.savefig(os.path.join(OUTDIR, f"{FIELD}_anomaly.png", dpi=300))
