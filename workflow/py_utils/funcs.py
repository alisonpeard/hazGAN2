"""
Functions to create variables of interest from ERA5 data.

Add custom functions as needed, inputs must be xarray DataArrays of
ERA5/IMDAA variables, outputs must be xarray DataArrays the same shape.

Reference here: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview

"""
import numpy as np
import xarray as xr

def unpack(params:xr.Dataset):
  """Unpack variables as named dict"""
  return {var: params[var] for var in params.data_vars}

def identity(x, *args, **kwargs):
    return x

def wind_speed(u, v, *args, **kwargs):
    return np.sqrt(u**2 + v**2)

def wind_direction(u, v, *args, **kwargs):
    return (np.arctan2(u, v) * 180 / np.pi + 360) % 360

# wind power functions
def wind_power(u10:xr.DataArray, v10:xr.DataArray,
               params:xr.Dataset) -> xr.DataArray:
  """
  Calculate wind power potential.

  Args:
    u10 : xr.DataArray
        10m u component of wind
    v10 : xr.DataArray
        10m v component of wind
  	params : xr.Dataset
    	Dataset with variable names corresponding to parameters for functions.
        For this function need:
        	- hub scaling parameters (hub_height, roughness)
            - 3PLE parameters (alpha, beta, ws0)
  
  References:
  ..[1] https://doi.org/10.1016/j.rser.2019.04.065
  ..[2] https://doi.org/10.1016/j.epsr.2017.10.028
  """
  def rescale_winds(ws, hub_height, gamma, h0=10.0, **kwargs):
    ws_rescaled = ws * (hub_height / h0)**gamma
    return ws_rescaled
  
  def ple3(ws, ws0, alpha, beta, **kwargs):
    power = alpha / (1 + np.exp(-beta * (ws - ws0)))
    return power
  
  ws = wind_speed(u10, v10) # should already be defined
  param_dict = unpack(params)
  ws_at_hub = rescale_winds(ws, **param_dict)
  potential = ple3(ws_at_hub, **param_dict)
  return potential


def solar_power(t2m:xr.DataArray, ssrd:xr.DataArray, u10:xr.DataArray,
                v10:xr.DataArray, params:xr.Dataset) -> xr.DataArray:
  """Calculate solar power potential.
  
  Args:
    t2m : xr.DataArray
        2m temperature
    ssrd : xr.DataArray
        Surface solar radiation downwards
    u10 : xr.DataArray
        10m u component of wind
    v10 : xr.DataArray
        10m v component of wind
    params : xr.Dataset
        Dataset with variable names corresponding to parameters for functions.
        For this function need:
            - solar cell temperature parameters (c_0, c_1, c_2, c_3)
            - performance ratio parameters (g_stc, t_ref, gamma)
    
    References:
    ..[1] https://doi.org/10.1016/j.rser.2019.04.065
    """
  def solar_cell_temp(t2m, sr, ws, c_0, c_1, c_2, c_3, **kwargs):
    return c_0 + c_1*t2m + c_2*sr + c_3*ws
  
  def performance_ratio(cell_temp, ref_temp, gamma, **kwargs):
    return 1 + gamma * (cell_temp - ref_temp)

  ws = wind_speed(u10, v10)
  param_dict = unpack(params)
  cell_temp = solar_cell_temp(t2m, ssrd, ws, **param_dict)
  pr = performance_ratio(cell_temp, **param_dict)
  potential = pr * (ssrd / params['g_stc'])
  return potential







