# Step 1: Mangrove damage probability fields
# convert samples to netcdf and add dependence assumption fields
# load samples as netcdf (including dependentce assumptions)
# add monthly medians
# load traning data -> 
# get yearly rate
# define damage ufunc and apply it to each netcdf samples, train, test
# NEW: save damages

# Step 2: Intersect damage fields with mangrove fields
# load mangroves
# intersect mangroves with damage fields
# get mangrove damages
# Calculate damagearea

# %%
import numpy as np
import os
import xarray as xr
from environs import Env
import matplotlib.pyplot as plt

from hazGAN.mangrove_demo import mangroveDamageModel
from analysis import load_samples, get_monthly_medians, λ

env = Env()
env.read_env()

# %% load generated data
THRESHOLD = 15. # None for all storms
TYPE = "trunc-1_0"
MODEL     = 24
MODEL     = str(MODEL).zfill(5) if isinstance(MODEL, int) else MODEL

samples_dir = env.str("SAMPLES_DIR")
data_dir    = env.str("DATA_DIR")

samples = load_samples(samples_dir, data_dir, MODEL, threshold=THRESHOLD, sampletype=TYPE)
medians = get_monthly_medians(data_dir, "September")
data    = xr.open_dataset(os.path.join(data_dir, "data.nc"))
nobs    = data.sizes['time']
nyears = len(np.unique(data['time.year']))

# get event sets
train       = samples['training']['data'] + medians
valid       = samples['valid']['data'] + medians
fake        = samples['samples']['data'] + medians
independent = samples['assumptions']['independent'] + medians
dependent   = samples['assumptions']['dependent'] + medians
dependent_uniform = samples['assumptions']['dependent_u']
dependent_rp      = samples['assumptions']['dependent_rp']

#  check orientation
fig, ax = plt.subplots()
ax.imshow(dependent[0, :, :, 0])
ax.set_title("I should be upside down!")

# %% Turn them all into netcdf
ref_data = xr.open_dataset(os.path.join(data_dir, "data.nc"))
coords = ref_data.coords
lat = coords['lat'].data
lon = coords['lon'].data

def to_xarray(array):
    da = xr.DataArray(array, coords=[
    ('sample', range(len(array))),
    ('lat', lat),
    ('lon', lon),
    ('field', ['u10', 'tp', 'mslp'])
    ])
    return da

trainda = to_xarray(train)
validda = to_xarray(valid)
fakeda = to_xarray(fake)
independentda = to_xarray(independent)
dependentda = to_xarray(dependent)
dependent_uniformda = to_xarray(dependent_uniform)
dependent_rpda = to_xarray(dependent_rp)

# %% make a dataset using their names
train = xr.Dataset({'train': trainda})
valid = xr.Dataset({'valid': validda})
fake = xr.Dataset({'fake': fakeda})
independent = xr.Dataset({'independent': independentda})

# %% dependent data needs some extra steps
dependent   = xr.Dataset({'dependent': dependentda})
dependent['uniform'] = dependent_uniformda
dependent['return_period'] = dependent_rpda
# %%
dependent['uniform'].plot.hist()

# convert return periods to 1-d array
assert np.isclose(dependent['uniform'].std(dim=['lat', 'lon']).max(), 0)
assert np.isclose(dependent['return_period'].std(dim=['lat', 'lon']).max(), 0)
dependent['uniform'] = dependent['uniform'].mean(dim=['lat', 'lon'])
dependent['return_period'] = dependent['return_period'].mean(dim=['lat', 'lon', 'field'])

# %% Let's try RP maps here"
if False:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from hazGAN.statistics import ecdf


    def calculate_hazard_map(rp, da:xr.DataArray, yearly_rate:int) -> xr.Dataset:
        """https://docs.xarray.dev/en/stable/examples/apply_ufunc_vectorize_1d.html"""

        def rp_ufunc(x):
            sorting = np.argsort(x)
            x_sorted = x[sorting]
            survival_p = 1 - ecdf(x_sorted)(x_sorted)
            return_periods = 1 / (yearly_rate * survival_p)
            out = np.interp(rp, return_periods, x_sorted) # consider better interpolator in hazGAN.statistics.empirical
            return out

        da = da.copy()

        RPs = xr.apply_ufunc(
            rp_ufunc,
            da,
            input_core_dims=[['sample']],
            output_core_dims=[[]],
            exclude_dims=set(('sample',)), # dimensions allowed to change size, must be set!
            vectorize=True,                # loop over non-core dimensions,
            dask="parallelized",
            output_dtypes=[float]
            )
        
        return RPs

    rp_maps = []
    rps = [2, 20, 100, 250, 500, 1000] #, 500, 1000, 2000]
    for rp in rps:
        rp_map = calculate_hazard_map(rp, fake['fake'].sel(field='u10'), rate)
        rp_maps.append(rp_map)

    # merge them
    rp_maps = xr.concat(rp_maps, dim='rp')
    rp_maps = rp_maps.assign_coords({'rp': rps})

    # 
    cmap = plt.get_cmap("YlOrRd")
    cmap.set_over('black')
    cmap.set_under('white')
    levels = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55]

    fig, axs = plt.subplots(1, len(rps), figsize=(24, 4), subplot_kw={'projection': ccrs.PlateCarree()})

    for ax in axs:
    # cartopy features
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, edgecolor='black')
        ax.add_feature(cfeature.RIVERS, edgecolor='lightblue')

    # add axis below for colorbar
    fig.subplots_adjust(bottom=0.2)

    # add new axis for colorbar
    cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])
    norm = plt.cm.colors.BoundaryNorm(levels, cmap.N)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax,
                        orientation='horizontal', extend='both', aspect=0.1)
    cbar.set_label('Wind speed (m/s)')
    cbar.set_ticks(levels)
    cbar.set_ticklabels(levels)

    # plot hazard map
    for i, rp in enumerate(rps):
        rp_maps.sel(rp=rp).plot.contourf(cmap=cmap, add_colorbar=False, extend='both', levels=levels, ax=axs[i])
        axs[i].set_title(f"1-in-{rp}")
    fig.suptitle("Return period maps for train winds")

# %% check for nans
for ds in [train, valid, fake, independent, dependent]:
    assert ds.isnull().sum() == 0, "Nans found in dataset"
    # NaNs found here! Why?
# %% predict mangrove damages
model = mangroveDamageModel()
train_damages = model.predict(train, ["train"])
valid_damages = model.predict(valid, ["valid"])
fake_damages = model.predict(fake, ["fake"])
independent_damages = model.predict(independent, ["independent"])
dependent_damages = model.predict(dependent, ["dependent"])

# TODO: make sure dependent has the extra fields
# %% rename the <>_damages to damage_prob
train_damages = train_damages.rename({"train_damage": "damage_prob"})
valid_damages = valid_damages.rename({"valid_damage": "damage_prob"})
fake_damages = fake_damages.rename({"fake_damage": "damage_prob"})
independent_damages = independent_damages.rename({"independent_damage": "damage_prob"})
dependent_damages = dependent_damages.rename({"dependent_damage": "damage_prob"})

# %% make return period maps from fake winds
train_damages.to_netcdf(os.path.join(samples_dir, "mangrove_demo", "train_damages.nc"))
valid_damages.to_netcdf(os.path.join(samples_dir, "mangrove_demo", "valid_damages.nc"))
fake_damages.to_netcdf(os.path.join(samples_dir, "mangrove_demo", "fake_damages.nc"))
independent_damages.to_netcdf(os.path.join(samples_dir, "mangrove_demo", "independent_damages.nc"))
dependent_damages.to_netcdf(os.path.join(samples_dir, "mangrove_demo", "dependent_damages.nc"))

# %%