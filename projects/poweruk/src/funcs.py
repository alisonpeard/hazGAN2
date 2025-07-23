import numpy as np
import xarray as xr


def convert_360_to_180(ds):
    """Convert longitude from 0-360 to -180 to 180."""
    ds = ds.assign_coords(longitude=(ds.longitude + 180) % 360 - 180)
    return ds.sortby('longitude')


def direction(ds:xr.Dataset, u:str, v:str,
              params={}
              ) -> xr.DataArray:
    """
    Calculate meteorological wind direction from ERA5 u and v components.
    
    Returns direction wind is coming FROM in degrees:
    - 0° = North (wind from north)
    - 90° = East (wind from east)  
    - 180° = South (wind from south)
    - 270° = West (wind from west)
    """
    # Convert to meteorological wind direction (where wind comes FROM) 
    u = ds[u]
    v = ds[v]
    direction = (270 - np.arctan2(v, u) * 180 / np.pi) % 360
    return direction


def sum_30_days(ds:xr.Dataset, arg:str, params={}) -> xr.DataArray:
    """
    Sum an hourly variable over the last 30 days.
    """
    return ds[arg].rolling(time=720).sum().dropna(dim="time", how="all")