"""
Functions for pre-processing input data and calculating derived variables.
"""
import numpy as np
import xarray as xr


def wind_speed(ds:xr.Dataset, u_arg:str, v_arg:str,
               params={}, dim="time") -> xr.DataArray:
    u = ds[u_arg]
    v = ds[v_arg]
    return np.sqrt(u**2 + v**2).max(dim=dim)


def wind_direction(u:xr.DataArray, v:xr.DataArray, *args, **kwargs) -> xr.DataArray:
    return (np.arctan2(u, v) * 180 / np.pi + 360) % 360