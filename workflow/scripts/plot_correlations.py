# %% 
import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from warnings import warn

import logging

from hazGAN.plotting import fields
from hazGAN.plotting import spatial

# set font to Helvetica
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Helvetica'

if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
        )
    TRAIN = snakemake.input.train
    GENER = snakemake.input.generated
    DIR0  = snakemake.output.dir0
    DIR1  = snakemake.output.dir1
    DATASET = snakemake.params.dataset.upper()
    DO_SUBSET = snakemake.params.do_subset
    THRESH = snakemake.params.event_subset
    OUTRES = snakemake.params.outres

    # TRAIN = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/training/data.nc"
    # GENER = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/generated/data.nc"
    # DIR0 = "/Users/alison/Local/hazGAN2/results/figures/correlations_field/"
    # DIR1 = "/Users/alison/Local/hazGAN2/results/figures/correlations_spatial/"
    # DATASET = "imdaa".upper()
    # DO_SUBSET   = True
    # THRESH = {
    #     "field": "ws",
    #     "func": "max",
    #     "value": 25.5,
    # }

    os.makedirs(DIR0, exist_ok=True)
    os.makedirs(DIR1, exist_ok=True)

    # load data and samples
    train = xr.open_dataset(TRAIN)
    if DO_SUBSET:
        train['intensity'] = getattr(train.sel(field=THRESH["field"]).anomaly, THRESH["func"])(dim=['lon', 'lat'])
        mask = (train['intensity'] > THRESH["value"]).values
        idx  = np.where(mask)[0]
        train   = train.isel(time=idx)
    
    train_x = train["anomaly"].values
    train_u = train["uniform"].values

    gener = xr.open_dataset(GENER)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values

    for FIELDS in [[0, 1], [0, 2], [1, 2]]:
        logging.info(f"Plotting fields {FIELDS[0]} and {FIELDS[1]}")
        # plot Smith (1990) extremal dependence
        fig_smith = fields.plot(gener_u, train_u, fields.smith1990, fields=FIELDS, figsize=.6,
                    title="", cbar_label=r"$\hat\theta$",
                    cmap="Spectral", vmin=1, vmax=4)
        fig_smith.savefig(os.path.join(DIR0, f"smith1990_{FIELDS[0]}-{FIELDS[1]}.png"), dpi=300)

        # plot Pearson correlation
        fig_pears = fields.plot(gener_u, train_u, fields.pearson, fields=FIELDS, figsize=.6,
                    title="", cbar_label=r"$r$", vmin=-1, vmax=1, cmap="Spectral_r")

        fig_pears.savefig(os.path.join(DIR0, f"pearson_{FIELDS[0]}-{FIELDS[1]}.png"), dpi=300)
    
    # %% - - - - - Plot spatial correlations - - - - - - - - - - - - - - - - -
    logging.info("Resampling to 16 x 16")
    outres = OUTRES
    ntrain, inx, iny, _ = train_u.shape
    ngener, _, _, _ = gener_u.shape
    bin_size = inx // outres
    train_u = train_u.reshape(ntrain, outres, bin_size, outres, bin_size, 3).max(3).max(1)
    gener_u = gener_u.reshape(ngener, outres, bin_size, outres, bin_size, 3).max(3).max(1)
    logging.info(f"Train shape: {train_u.shape}")
    logging.info(f"Generated shape: {gener_u.shape}")

    for FIELD in [0, 1, 2]:
        logging.info(f"Plotting field {FIELD}")
        
        fig_smith = spatial.plot(gener_u, train_u, spatial.smith1990, field=FIELD, figsize=.6,
                    title="", cbar_label=r"$\hat\theta$",
                    cmap="Spectral", vmin=1, vmax=4)
        fig_smith.savefig(os.path.join(DIR1, f"smith1990_{FIELD}.png"), dpi=300)
        # plot Pearson correlation
        fig_pears = spatial.plot(gener_u, train_u, spatial.pearson, field=FIELD, figsize=.6,
                    title="", cbar_label=r"$r$", vmin=-1, vmax=1, cmap="Spectral_r")
        fig_pears.savefig(os.path.join(DIR1, f"pearson_{FIELD}.png"), dpi=300)
    # %%