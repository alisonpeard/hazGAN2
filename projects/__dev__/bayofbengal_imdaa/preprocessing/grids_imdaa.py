# %%
import os
from glob import glob
import xarray as xr

# %%
imdaapath = os.path.join("/soge-home", "data", "analysis", "imdaa", "gust", "*.nc")
imdaafiles = glob(imdaapath)
ds = xr.open_dataset(imdaafiles[0])
ds = ds.isel(time=0)
print(ds)

ds = ds.rename({"GUST_10m": "sample_param"})
ds["sample_param"] = 1 + (ds["sample_param"] * 0)
compression = {'zlib': True, 'complevel': 5}
encoding = {var: compression for var in ds.data_vars}
ds.to_netcdf("../../resources/grids/imdaa.nc", engine='netcdf4', encoding=encoding)
ds.to_netcdf("../../resources/params/bayofbengal_imdaa.nc", engine='netcdf4', encoding=encoding)
print("Saved imdaa.nc grid file.")
# %%