# %% 
import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from warnings import warn

import logging

from hazGAN.plotting import scatter
from hazGAN.constants import OBSERVATION_POINTS
from hazGAN import op2idx # in utils.py, need to check method

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
    DIR  = snakemake.output.outdir
    DO_SUBSET = snakemake.params.do_subset
    THRESH = snakemake.params.event_subset
    YMIN = snakemake.params.ymin
    YMAX = snakemake.params.ymax
    XMIN = snakemake.params.xmin
    XMAX = snakemake.params.xmax
    POIS = snakemake.params.pois

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
    
    # 
    ops = op2idx(OBSERVATION_POINTS, train_x[0, ..., 0], extent=[80, 95, 10, 25])
    
    # TODO (optional): add Brown–Resnick data
    # brownresnick = pd.read_parquet(os.path.join(br_dir, "ECs_u10.parquet"))
    
    pixels = [ops['chittagong'], ops['dhaka']]

    for FIELD in range(3):
        fig = scatter.plot(gener_x, train_x, field=FIELD, pixels=pixels, s=10,
                        cmap='viridis', xlabel="Chittagong", ylabel="Dhaka")

        fig.savefig(os.path.join(DIR, f"scatterplot_{FIELD}.png"), dpi=300)
        plt.close(fig)
