# %%
from pathlib import Path
import glob as glob
import xarray as xr

input_data = "/Users/alison/Documents/dphil/data/hazGAN2/projects/poweruk_winter/results/training/data.nc"

ds = xr.open_dataset(input_data)
ds.uniform.shape
nimgs = ds.time.size
u = ds.uniform.values
# %%
print(u.max() < 1.)
print(u.min() > 0.)
# %%

import sys
sys.path.append("workflow")
import src.python.statistics as stats
import matplotlib.pyplot as plt

domain = "laplace"
ppf = getattr(stats, domain)

# %%
y = ppf(u)
# %%

plt.hist(y.flatten(), bins=100)
# %%
rp_max = 1e6
# %%
ymin = ppf(1 / rp_max)
ymax = ppf(1 - 1 / rp_max)
print(ymin, ymax)
# %%
y_scaled = (y - ymin) / (ymax - ymin)

# %%
