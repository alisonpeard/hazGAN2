# %%
import xarray as xr

path = "/Users/alison/Documents/dphil/data/hazGAN2/projects/poweruk_winter/results/processing/data.nc"

ds = xr.open_dataset(path)

u10 = ds.sel(field="u10_gust")
v10 = ds.sel(field="v10_gust")
r30 = ds.sel(field="r30")

# %%
freq = ds.attrs['yearly_freq']

for field in [u10, v10, r30]:
    pmax = field['uniform'].max().values
    print(f"{pmax:.8f} m/s", end=" - ")
    print(f"1-in-{1/(1-pmax):.0f} storms", end=" - ")
    print(f"1-in-{1/(freq*(1-pmax)):.0f} years")
# %%

# %%
shapemin = v10.sel(param="shape_upper")['params'].max()
shapemin

# %%