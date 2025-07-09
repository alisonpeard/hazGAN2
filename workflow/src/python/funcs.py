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


# the following functions are for deriving variables from the input data:
def identity(x:xr.Dataset, arg, **kwargs) -> xr.DataArray:
    return x[arg]


def max(ds:xr.Dataset, arg, **kwargs) -> xr.DataArray:
    """
    Calculate the maximum value of a variable.
    """
    return ds[arg].max(dim="time", skipna=True)


def arg2max(ds:xr.Dataset, arg1, arg2, **kwargs) -> xr.DataArray:
    """
    Calculate the index of the maximum value of a variable.
    """
    idx_max = ds[arg2].argmax(dim="time")
    return ds.isel(time=idx_max)[arg1]


def mean(ds:xr.Dataset, arg, **kwargs) -> xr.DataArray:
    """
    Calculate the mean value of a variable.
    """
    return ds[arg].mean(dim="time", skipna=True)