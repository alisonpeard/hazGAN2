# %% histograms
import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import logging
from warnings import warn

from hazGAN.plotting import misc

# set font to Helvetica
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Helvetica'

if __name__ == "__main__":
    TRAIN = snakemake.input.train
    GENER = snakemake.input.generated
    MONTH   = snakemake.params.month
    FIGURE  = snakemake.output.figure
    DATASET = snakemake.params.dataset.upper()
    DO_SUBSET = snakemake.params.do_subset
    THRESH = snakemake.params.event_subset

    # configure logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info(f"Plotting barcharts for {DATASET} dataset, month {MONTH}")

    # TRAIN = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/training/data.nc"
    # GENER = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/generated/data.nc"
    # FIGURE = "/Users/alison/Local/hazGAN2/results/figures/samples/barchart.png"
    # DATASET = "imdaa".upper()
    # DO_SUBSET   = True
    # THRESH = {
    #     "field": "ws",
    #     "func": "max",
    #     "value": 25.5,
    # }

    # load data and samples
    train = xr.open_dataset(TRAIN)
    medians = train.sel(month=[MONTH]).medians.values
    
    if DO_SUBSET:
        train['intensity'] = getattr(train.sel(field=THRESH["field"]).anomaly, THRESH["func"])(dim=['lon', 'lat'])
        mask = (train['intensity'] > THRESH["value"]).values
        idx  = np.where(mask)[0]
        train   = train.isel(time=idx)
    
    train_x = train["anomaly"].values

    gener = xr.open_dataset(GENER)
    gener_x = gener["anomaly"].values

    # add medians for month back to samples
    logging.info(f"{medians.shape=}, {train_x.shape=}, {gener_x.shape=}")
    train_x += medians
    gener_x += medians

    bar_width = 0.25
    fig, ax = misc.saffirsimpson_barchart(gener_x, train_x, bar_width=bar_width, title="")

    if False: # load and add IBTrACS
        raise NotImplementedError("IBTrACS not implemented")
        from hazGAN import saffirsimpson

        ibtracs = pd.read_csv("/Users/alison/Documents/DPhil/paper1.nosync/training/18x22/ibtracs_dates.csv")
        ibtracs["time"] = pd.to_datetime(ibtracs["time"])
        ibtracs = ibtracs.sort_values("time", ascending=False)
        ibtracs = ibtracs[ibtracs['event'] != "Not_Named"]
        ibtracs_storms = ibtracs.groupby("event").agg({"time": ["min", "max"], "wind": "max"})
        ibtracs_storms.columns = ["start", "end", "wind"]
        stormibtracs_stormss = ibtracs_storms.sort_values("start", ascending=False)
        ibtracs_storms = ibtracs_storms.sort_values("wind", ascending=False)
        nstorms = len(ibtracs_storms)
        ibtracs_storms['Category'] = ibtracs_storms['wind'].apply(saffirsimpson)
        cat_counts = ibtracs_storms['Category'].value_counts().sort_index()
        cat_density = cat_counts / len(ibtracs_storms)
        all_categories = np.arange(-1, 6)
        ibtracs_density = pd.Series(cat_density.get(cat, 0) for cat in all_categories)
        r = np.arange(len(all_categories))
        r = [x + 2 * bar_width + 0.01 for x in r]

        ax.bar(r, ibtracs_density, width=bar_width, color="#F6F5EE", label="IBTrACS",
            edgecolor='#666666', linewidth=0.5)
        
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=300)
# %%
