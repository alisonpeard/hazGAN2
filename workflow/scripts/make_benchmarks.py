
"""Create total dependence/independence benchmarks."""
# %%
import os
from glob import glob
import logging

os.environ["USE_PYGEOS"] = "0"
from PIL import Image
import numpy as np
import xarray as xr

from hazGAN.statistics import gumbel, inv_gumbel, invPIT


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
    
    # load parameters
    INPUT   = snakemake.input.image_dir
    RESX    = snakemake.params.resx
    RESY    = snakemake.params.resy
    YEAR0   = snakemake.params.year0
    YEARN   = snakemake.params.year1
    DO_SUBSET = snakemake.params.do_subset
    THRESH = snakemake.params.event_subset
    NYRS    = snakemake.params.nyrs # number of years of data to generate
    FIELDS = snakemake.params.fields
    data = xr.open_dataset(snakemake.input.data)
    N_HAZMAPS = 10

    # %%
    params = data.params.values
    distns = [field["distn"] for field in FIELDS]
    ref = data.copy()
    if DO_SUBSET:
        data['intensity'] = getattr(data.sel(field=THRESH["field"]).anomaly, THRESH["func"])(dim=['lon', 'lat'])
        mask = (data['intensity'] > THRESH["value"]).values
        idx  = np.where(mask)[0]
        data   = data.isel(time=idx)

    nyears = YEARN - YEAR0
    h, w, c = RESY, RESX, 3

    x_ref = ref.anomaly.values
    x = data.anomaly.values
    # %%
    λ_event  = λ(len(x_ref), nyears) # event rate for all storms
    λ_extreme = λ(len(x), nyears) # event rate for extreme storms
    nevents = int(NYRS * λ_event)  # how many storms for NYRS years of data
    nextremes = int(NYRS * λ_extreme) # how many extreme storms for NYRS years of data

    # %%
    # make totally independent samples
    # sample independent events from base distribution
    independent_u = np.random.uniform(1e-6, 1-1e-6, size=(nevents, h, w, c))

    # make hazard map-style (dependent) samples
    dependent_rps = np.logspace(0, 3, num=N_HAZMAPS, base=10)
    dependent_survival  = 1 / (λ_event * dependent_rps)
    dependent_u   = 1 - dependent_survival
    dependent_u   = np.repeat(dependent_u, h*w*c, axis=0)
    dependent_u   = dependent_u.reshape(N_HAZMAPS, h, w, c)
    dependent_rps       = np.repeat(dependent_rps, h*w*c, axis=0)
    dependent_rps       = dependent_rps.reshape(N_HAZMAPS, h, w, c)

    # %% transform to data space

    independent_x = invPIT(independent_u, x_ref, params, distns=distns)
    dependent_x   = invPIT(dependent_u, x_ref, params, distns=distns)

