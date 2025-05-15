
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
    """Calculate the extreme event rate for a given return period."""
    return number_extremes / (total_years / unit)

if __name__ == "__main__":
    # configure logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    
    # load i/o
    DATA   = snakemake.input.data
    DEPENDENT = snakemake.output.dependent
    INDEPENDENT = snakemake.output.independent

    # load parameters
    RESX    = snakemake.params.resx
    RESY    = snakemake.params.resy
    YEAR0   = snakemake.params.year0
    YEARN   = snakemake.params.yearn
    NYRS    = snakemake.params.nyrs # number of years of data to generate
    FIELDS = snakemake.params.fields
    N_HAZMAPS = snakemake.params.n_hazmaps

    # load the training data
    data = xr.open_dataset(DATA)
    params = data.params.values
    distns = [field["distn"] for field in FIELDS]

    nyears = YEARN - YEAR0
    h, w, c = RESY, RESX, 3

    x = data.anomaly.values

    # calculate event rate and number of events to generatee
    λ_event  = λ(len(x), nyears)   # event rate for all storms
    nevents = int(NYRS * λ_event)  # how many storms for NYRS years of data

    # make totally independent samples (sample from base distribution)
    independent_u = np.random.uniform(1e-6, 1-1e-6, size=(nevents, h, w, c))

    # make hazard map-style (dependent) samples
    dependent_rps = np.logspace(0, 3, num=N_HAZMAPS, base=10)
    dependent_survival  = 1 / (λ_event * dependent_rps)
    dependent_u   = 1 - dependent_survival
    dependent_u   = np.repeat(dependent_u, h*w*c, axis=0)
    dependent_u   = dependent_u.reshape(N_HAZMAPS, h, w, c)
    dependent_rps       = np.repeat(dependent_rps, h*w*c, axis=0)
    dependent_rps       = dependent_rps.reshape(N_HAZMAPS, h, w, c)

    # transform to data space
    independent_x = invPIT(independent_u, x, params, distns=distns)
    dependent_x   = invPIT(dependent_u, x, params, distns=distns)

    # make datasets
    independent_ds = xr.Dataset(
        data_vars={
            "anomaly": (("time", "lat", "lon"), independent_x),
            "uniform": (("time", "lat", "lon"), independent_u),
        },
        coords={
            "time": np.arange(nevents),
            "lat": data.lat,
            "lon": data.lon,
        },
    )
    independent_ds["params"] = (("time", "param"), params)

    dependent_ds = xr.Dataset(
        data_vars={
            "anomaly": (("time", "lat", "lon"), dependent_x),
            "uniform": (("time", "lat", "lon"), dependent_u),
        },
        coords={
            "time": np.arange(nevents),
            "lat": data.lat,
            "lon": data.lon,
        },
    )
    dependent_ds["params"] = (("time", "param"), params)
    dependent_ds["rps"] = (("time", "hazmap"), dependent_rps)

    logging.info(f"Dependent dataset shape: {dependent_ds.anomaly.shape}")
    logging.info(f"Independent dataset shape: {independent_ds.anomaly.shape}")

    # save datasets
    independent_ds.to_netcdf(INDEPENDENT, mode="w", format="NETCDF4")
    dependent_ds.to_netcdf(DEPENDENT, mode="w", format="NETCDF4")
    logging.info(f"Saved independent dataset to {INDEPENDENT}")
    logging.info(f"Saved dependent dataset to {DEPENDENT}")

