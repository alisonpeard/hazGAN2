# %%
import os
import numpy as np
import xarray as xr

u10file = "/Users/alison/Documents/DPhil/data/era5/10m_u_component_of_wind/nc/10m_u_component_of_wind_2020.nc"
v10file = "/Users/alison/Documents/DPhil/data/era5/10m_v_component_of_wind/nc/10m_v_component_of_wind_2020.nc"
gustfile = "/Users/alison/Documents/DPhil/data/era5/instantaneous_10m_wind_gust/nc/instantaneous_10m_wind_gust_2020.nc"
tpfile = "/Users/alison/Documents/DPhil/data/era5/total_precipitation/nc/total_precipitation_2020.nc"

gust = xr.open_dataset(gustfile, engine='netcdf4')
u10 = xr.open_dataset(u10file, engine='netcdf4')
v10 = xr.open_dataset(v10file, engine='netcdf4')
tp = xr.open_dataset(tpfile, engine='netcdf4')

# rename valid_time to time
u10 = u10.rename({"valid_time": "time"})
v10 = v10.rename({"valid_time": "time"})
gust = gust.rename({"valid_time": "time"})
tp = tp.rename({"valid_time": "time"})
# %%
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

# %%
temp = relative_direction(u10.u10, v10.v10, gust.i10fg)
res = temp.resample(time='1D').map(dir_at_max_wind)

import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(15, 6))
temp.isel(dir_speed=0, time=0).plot(cmap='twilight', ax=axs[0])
temp.isel(dir_speed=1, time=0).plot(cmap='viridis', vmin=0, vmax=30, ax=axs[1])
res.isel(time=0).plot(cmap='twilight', ax=axs[2])
# %%
tp_rolling = sum_n_hours(tp.tp)
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
tp.tp.isel(time=0).plot(cmap='PuBu', ax=axs[0])
tp_rolling.isel(time=0).plot(cmap='PuBu', ax=axs[1])

# %%
