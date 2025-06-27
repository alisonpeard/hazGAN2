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


def identity(x:xr.DataArray, *args, **kwargs) -> xr.DataArray:
    return x