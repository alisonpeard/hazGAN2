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
    - 0° = northerly (wind from north)
    - 90° = easterly (wind from east)  
    - 180° = southerly (wind from south)
    - 270° = westerly (wind from west)

    See:
    - https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
    - https://numpy.org/doc/2.2/reference/generated/numpy.arctan2.html
    """
    # Convert to meteorological wind direction (where wind comes FROM)
    u = ds[u]
    v = ds[v]
    direction = (180 + 180 / np.pi * np.arctan2(v, u)) % 360
    return direction


def sum_30_days(ds:xr.Dataset, arg:str, params={}) -> xr.DataArray:
    """
    Sum an hourly variable over the last 30 days. Remove NaNs.
    """
    t0 = ds["time"].min().data
    tn = ds["time"].max().data
    t30 = t0 + np.timedelta64(30, 'D')
    ds = ds.sortby("time")

    r30 = ds[arg].rolling(time=720).sum()
    r30 = r30.sel(time=slice(t30, tn))

    return r30


def dir_at_max_gust(
    ds:xr.Dataset, u:str, v:str, i10fg:str,
    params={}
    ) -> xr.Dataset:
    """Calculate wind direction at time of maximum gust."""
    gust = ds[i10fg]
    u    = ds[u]
    v    = ds[v]
    tmax = gust.idxmax(dim="time")
    u_max = u.sel(time=tmax)
    v_max = v.sel(time=tmax)
    direction_max = (180 + 180 / np.pi * np.arctan2(v_max, u_max)) % 360
    return direction_max


def scale_to_gust(
    ds:xr.Dataset, u:str, v:str, i10fg:str,
    params={}
    ) -> xr.Dataset:
    """Scale the u and v components of wind to gust speed.
    
    This way, wind speeds in analysis will be the ERA5 wind
    speeds, but wind direction will be unchanged.
    """
    gust = ds[i10fg]
    u    = ds[u]
    v    = ds[v]
    gust_factor = gust / np.sqrt(u**2 + v**2)
    u_gust = gust_factor * u
    return u_gust
