import numpy as np
import xarray as xr



def direction(u: xr.DataArray, v: xr.DataArray,
              *args, **kwargs) -> xr.DataArray:
    """
    Calculate meteorological wind direction from ERA5 u and v components.
    
    Returns direction wind is coming FROM in degrees:
    - 0° = North (wind from north)
    - 90° = East (wind from east)  
    - 180° = South (wind from south)
    - 270° = West (wind from west)
    """
    # Convert to meteorological wind direction (where wind comes FROM)
    direction = (270 - np.arctan2(v, u) * 180 / np.pi) % 360
    return direction


def sum_n_hours(tp:xr.DataArray, n=720, *args, **kwargs) -> xr.DataArray:
    """
    Sum precipitation over n hours.
    """
    if tp.time.size < n:
        raise ValueError(f"Time dimension size {tp.time.size} is less than n={n}.")
    return tp.rolling(time=n).sum().dropna(dim="time", how="all")


def max(ds_grouped:xr.Dataset, *args) -> xr.DataArray:
    """
    Calculate the maximum value of a variable.
    """
    return ds_grouped[args[0]].max(dim="time", skipna=True)


def arg2max(ds_grouped:xr.Dataset, *args) -> xr.DataArray:
    """
    Calculate the index of the maximum value of a variable.
    """
    idx_max = ds_grouped[args[1]].argmax(dim="time")
    return ds_grouped.isel(time=idx_max)[args[0]]


def mean(ds_grouped:xr.Dataset, *args) -> xr.DataArray:
    """
    Calculate the mean value of a variable.
    """
    return ds_grouped[args[0]].mean(dim="time", skipna=True)
