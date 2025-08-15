# %%
import os
import numpy as np
import xarray as xr

wd = os.path.join("..", "..", "results")
train_path = os.path.join(wd, "training", "data.nc")

"""Lets find a number to extract 100 from the ~400 samples here."""
ds = xr.open_dataset(train_path)

def l2norm_max(ds:xr.Dataset, arg1:str, arg2:str,
              params={}, dim="time") -> xr.DataArray:
    l2norm = np.sqrt(ds[arg1]**2 + ds[arg2]**2)
    return l2norm.max(dim=dim)




func = l2norm_max
args = ["u10", "v10"]
for arg in args:
    ds[arg] = ds.sel(field=arg).anomaly

intensity = func(ds, *args, dim=["lon", "lat"])

for arg in args:
    ds = ds.drop_vars(arg)

# %%
import matplotlib.pyplot as plt
plt.hist(intensity.values.flatten(), bins=100);

# %%
sum((intensity > 25))
# %%
mask = (intensity > 25).values