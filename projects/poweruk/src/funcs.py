import numpy as np
import xarray as xr


def relative_direction(u:xr.DataArray, v:xr.DataArray,
                       gust:xr.DataArray, *args, **kwargs) -> xr.DataArray:
    """Deviation of wind direction from prevailing direction."""
    direction = (np.arctan2(u, v) * 180 / np.pi + 360) % 360
    prevailing_direction = direction.mean(dim=["time"], skipna=True)
    relative_direction = (direction - prevailing_direction) % 360
    temp = xr.concat([relative_direction, gust], dim="dir_speed")
    return temp


def dir_at_max_wind(temp:xr.DataArray, *args, **kwargs) -> xr.DataArray:
    """
    Calculate the direction at maximum gust speed.
    """
    direction = temp.sel(dir_speed=0)
    idx_gust = temp.isel(dir_speed=1).argmax(dim="time")
    direction = direction.isel(time=idx_gust)
    return direction


def sum_n_hours(tp:xr.DataArray, n=720, *args, **kwargs) -> xr.DataArray:
    """
    Sum precipitation over n hours.
    """
    if tp.time.size < n:
        raise ValueError(f"Time dimension size {tp.time.size} is less than n={n}.")
    return tp.rolling(time=n).sum().dropna(dim="time", how="all")