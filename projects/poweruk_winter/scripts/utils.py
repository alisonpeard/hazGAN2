import numpy as np
import xarray as xr


def handle_expver(ds):
    for coord in ["expver", "number"]:
        if coord in ds.dims:
            ds = ds.sel({coord: ds[coord][0]})
        if coord in ds.coords:
            ds = ds.reset_coords(coord, drop=True)
    return ds


def check_missing_days(ds):
    all_days = np.arange(ds.time.min().values, ds.time.max().values + np.timedelta64(1, 'D'), dtype='datetime64[D]')
    missing_days = np.setdiff1d(all_days, ds.time.values)
    if len(missing_days) > 0:
        print(f"Missing days: {len(missing_days)}")
        print(missing_days)
    else:
        print("No missing days found.")



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
    u = ds[u]
    v = ds[v]
    direction = (180 + 180 / np.pi * np.arctan2(v, u)) % 360
    return direction


def derive_variables(ds):
    """
    Derive variables from the dataset.
    """
    ds = handle_expver(ds)
    
    # Calculate wind speed and direction
    ds["vx"] = np.sqrt(ds.u10**2 + ds.v10**2)
    ds["dx"] = direction(ds, "u10", "v10")

    return ds