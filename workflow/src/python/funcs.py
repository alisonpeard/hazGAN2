"""
Functions for pre-processing input data and calculating derived variables.

Add custom functions as needed, inputs must be xarray DataArrays of
ERA5/IMDAA variables, outputs must be xarray DataArrays the same shape.

Reference here: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview

"""
import numpy as np
import xarray as xr


def unpack(params:xr.Dataset):
  """Unpack variables as named dict"""
  return {var: params[var] for var in params.data_vars}


# the following `init` functions are for deriving variables from the input data:
def identity(ds:xr.Dataset, arg:str, params={}) -> xr.DataArray:
    """Return the variable as is."""
    return ds[arg]


# these are hfuncs for aggregating across a time/sample dimension
def max(ds:xr.Dataset, arg:str, params={}, dim="time") -> xr.DataArray:
    """Calculate the maximum value of a variable."""
    return ds[arg].max(dim=dim, skipna=True)


def l2norm_max(ds:xr.Dataset, arg1:str, arg2:str,
              params={}, dim="time") -> xr.DataArray:
    l2norm = np.sqrt(ds[arg1]**2 + ds[arg2]**2)
    return l2norm.max(dim=dim)


def arg2max(ds:xr.Dataset, arg1:str, arg2:str, params={}, dim="time") -> xr.DataArray:
    """Calculate the index of the maximum value of a variable."""
    arg2_max = ds[arg2].max(dim=dim)
    mask = ds[arg2] == arg2_max
    return ds[arg1].where(mask).max(dim=dim)


def l2norm_argmax(ds:xr.Dataset, arg1:str, arg2:str,
              params={}, dim="time") -> xr.DataArray:
    l2norm = np.sqrt(ds[arg1]**2 + ds[arg2]**2)
    l2norm_max = l2norm.max(dim=dim)
    mask = l2norm == l2norm_max
    return ds[arg1].where(mask).max(dim=dim)


def mean(ds:xr.Dataset, arg, params={}, dim="time") -> xr.DataArray:
    """Calculate the mean value of a variable."""
    return ds[arg].mean(dim=dim, skipna=True)