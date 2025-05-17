
"""Create total dependence/independence benchmarks."""
# %%
import os
import logging

os.environ["USE_PYGEOS"] = "0"
from PIL import Image
import numpy as np
import xarray as xr

from hazGAN.statistics import invPIT


def λ(number_extremes:int, total_years:int, unit:int=1) -> float:
    """Calculate the event rate for a given return period unit (e.g., years)."""
    return number_extremes / (total_years / unit)

# %%
if __name__ == "__main__":
    # configure logging
    # logging.basicConfig(
    #     filename=snakemake.log.file,
    #     level=logging.INFO,
    #     format="%(asctime)s - %(levelname)s - %(message)s",
    # )
    
    # load i/o
    TRAIN   = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/training/data.nc"
    DEPENDENT = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/training/dependent.nc"
    INDEPENDENT = "/Users/alison/Local/hazGAN2/results/bayofbengal_imdaa/training/independent.nc"

    # load parameters
    RESX    = 64
    RESY    = 64
    YEAR0   = 1979
    YEARN   = 2023
    NYRS    = 1000 # number of years of data to generate
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

    # load the training data
    data   = xr.open_dataset(TRAIN)
    params = data["params"].values
    distns = [field["distn"] for field in FIELDS.values()]

    nyears  = YEARN - YEAR0
    h, w, c = RESY, RESX, 3

    x = data.anomaly.values

    # %% calculate event rate and number of events to generate
    NYRS = 500
    λ_event  = λ(len(x), nyears)   # event rate for all storms
    nevents = int(NYRS * λ_event)  # how many storms for NYRS years of data
    logging.info(f"Event rate: {λ_event:.2f} events/year")
    logging.info(f"Number of events: {nevents} events")

    logging.info(f"Generating {nevents} events for {NYRS} years of data")
    logging.info(f"Generating tensor of size {nevents} x {h} x {w} x {c} = {nevents * h * w * c:,.0f} total elements")

    # %% make totally independent samples (sample from base distribution)
    independent_u = np.random.uniform(1e-6, 1-1e-6, size=(nevents, h, w, c))

    # %%
    independent_x = invPIT(independent_u, x, params, distns=distns) # slow

    # %%
    independent_ds = xr.Dataset(
        data_vars={
            "anomaly": (("time", "lat", "lon", "field"), independent_x),
            "uniform": (("time", "lat", "lon", "field"), independent_u),
            "params": (("lat", "lon", "param", "field"), params),
        },
        coords={
            "field": list(FIELDS.keys()),
            "param": ["loc", "scale", "shape"],
            "time": np.arange(nevents),
            "lat": data.lat.values,
            "lon": data.lon.values,
        },
    )
    logging.info(f"Finished independent dataset")
    # %%
    independent_ds.isel(field=0, time=0).anomaly.plot()
    # %% make hazard map-style (dependent) samples
    dependent_rps = np.array([1, 2, 5, 10, 25, 50, 100, 250, 500, 1000])
    dependent_survival  = 1 / (λ_event * dependent_rps)
    dependent_u   = 1 - dependent_survival
    dependent_u   = np.repeat(dependent_u, h*w*c, axis=0)
    dependent_u   = dependent_u.reshape(len(dependent_rps), h, w, c)
    dependent_x   = invPIT(dependent_u, x, params, distns=distns)

    dependent_ds = xr.Dataset(
        data_vars={
            "anomaly": (("rp", "lat", "lon", "field"), dependent_x),
            "uniform": (("rp", "lat", "lon", "field"), dependent_u),
            "params": (("lat", "lon", "param", "field"), params)
        },
        coords={
            "field": list(FIELDS.keys()),
            "param": ["loc", "scale", "shape"],
            "rp": dependent_rps,
            "lat": data["lat"].values,
            "lon": data["lon"].values
        },
    )
    logging.info(f"Finished dependent dataset")
    # dependent_ds.isel(field=0, rp=9).anomaly.plot()

    # %% save datasets
    logging.info(f"Dependent dataset shape: {dependent_ds.anomaly.shape}")
    logging.info(f"Independent dataset shape: {independent_ds.anomaly.shape}")

    independent_ds.to_netcdf(INDEPENDENT)
    dependent_ds.to_netcdf(DEPENDENT)

    logging.info(f"Saved independent dataset to {INDEPENDENT}")
    logging.info(f"Saved dependent dataset to {DEPENDENT}")


# %%
