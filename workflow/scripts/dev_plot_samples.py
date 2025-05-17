"""Make 64x64 plots of data"""
# %%
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
# from snakemake.scripts import snakemake
from hazGAN.plotting import samples
from hazGAN.statistics import gumbel


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'

# %%
if __name__ == "__main__":
    TRAIN = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/training/data.nc"
    GENER = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/generated/data.nc"
    SHUFFLE = True
    OUTDIR  = "/Users/alison/Local/hazGAN2/results/figures/samples"
    DATASET = "imdaa"
    DO_SUBSET   = True
    THRESH = {
        "field": "ws",
        "func": "max",
        "value": 25.5,
    }

    FIELDS = {
        "ws": {
            "args": ["GUST_10m"],
            "func": "identity",
            "hfunc": "max",
            "obj": "max",
            "distn": "weibull",
            "units": "mps",
            "cmap": "viridis"
        },
        "msl": {
            "args": ["PRMSL_msl"],
            "func": "identity",
            "hfunc": "min",
            "obj": "min",
            "distn": "genpareto",
            "units": "Pa",
            "cmap": "Spectral"
        },
        "tp": {
            "args": ["APCP_sfc"],
            "func": "identity",
            "hfunc": "sum",
            "obj": "max",
            "distn": "genpareto",
            "units": "mm",
            "cmap": "PuBu"
        },
    }

    # load data and samples
    train = xr.open_dataset(TRAIN)
    if DO_SUBSET:
        mask = getattr(train.sel(field=THRESH["field"]).anomaly, THRESH["func"])(dim=["lon", "lat"]) > THRESH["value"]
        mask = np.where(mask.values)[0].tolist()
        train =train.isel(time=mask)
    
    train_x = train["anomaly"].values
    train_u = train["uniform"].values
    train_g = gumbel(train_u)

    gener = xr.open_dataset(GENER)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values
    gener_g = gumbel(gener_u)

    # %%
    if SHUFFLE:
        train_ids = np.random.permutation(train_u.shape[0])
        gener_ids = np.random.permutation(gener_u.shape[0])
    else:
        train_ids = np.arange(train_u.shape[0])
        gener_ids = np.arange(gener_u.shape[0])

    # %% make the samples
    for i, FIELD in enumerate(FIELDS.keys()):
        METRIC = FIELDS[FIELD]["units"]
        CMAP   = FIELDS[FIELD]["cmap"]

        figa = samples.plot(gener_g[gener_ids], train_g[train_ids], field=i, title="", cmap=CMAP, ndecimals=0)
        figb = samples.plot(gener_u[gener_ids], train_u[train_ids], field=i, title="", cbar_label="", cmap=CMAP, ndecimals=1)
        figc = samples.plot(gener_x[gener_ids], train_x[train_ids], field=i, title="", cbar_label=METRIC, cmap=CMAP)#, alpha=1e-6);

        os.makedirs(OUTDIR, exist_ok=True)

        figa.savefig(os.path.join(OUTDIR, f"{FIELD}_gumbel.png"), dpi=300)
        figb.savefig(os.path.join(OUTDIR, f"{FIELD}_uniform.png"), dpi=300)
        figc.savefig(os.path.join(OUTDIR, f"{FIELD}_anomaly.png"), dpi=300)

# %%
