"""
Functions for pre-processing input data and calculating derived variables.
"""
import numpy as np
import xarray as xr


def wind_speed(u:xr.DataArray, v:xr.DataArray, *args, **kwargs) -> xr.DataArray:
    return np.sqrt(u**2 + v**2)


def wind_direction(u:xr.DataArray, v:xr.DataArray, *args, **kwargs) -> xr.DataArray:
    return (np.arctan2(u, v) * 180 / np.pi + 360) % 360