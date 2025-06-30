"""
Step 2: Intersect damage fields with mangrove fields
load mangroves
intersect mangroves with damage fields
get mangrove damages
Calculate damagearea
"""
# %%
import os
import yaml
import xarray as xr
import matplotlib.pyplot as plt

from .mangroves.statistics import calculate_total_return_periods # calculate_eads

if __name__ == "__main__":
    # load project configuration
    with open(os.path.join("..", "config.yaml"), 'r') as stream:
        config = yaml.safe_load(stream)
    
    wd = os.path.join("..", "results", "mangroves")
    nyrs = config["yearn"] - config["year0"]

    # load damage data
    THRESHOLD = config['event_subset']['threshold']
    train_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "train.nc"))
    gener_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "gener.nc"))
    indep_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "indep.nc"))
    depen_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "depen.nc"))
    mangr_grid    = xr.open_dataset(os.path.join(wd, "mangrove_grid.nc"))

    # expected damage for each grid cell
    train_damages['expected_damage'] = train_damages['damage_prob'] * mangr_grid['area']
    gener_damages['expected_damage'] = gener_damages['damage_prob'] * mangr_grid['area']
    indep_damages['expected_damage'] = indep_damages['damage_prob'] * mangr_grid['area']
    depen_damages['expected_damage'] = depen_damages['damage_prob'] * mangr_grid['area']

    # make a datatree with return period variables
    def calculate_total_return_periods(
            damages:xr.Dataset, var:str
            ) -> xr.Dataset:
        if len(damages.data_vars) > 0: # skip root node in datatree
            # aggregate to overall damages
            npy = damages['rate'].values.item()
            totals = damages[var].sum(dim=['lat', 'lon']).to_dataset()
            
            # calculate return periods
            N = totals[var].sizes['sample']
            rank = totals[var].rank(dim='sample')
            totals['exceedence_prob'] = 1 - rank / (1 + N)
            totals['return_period'] = 1 / (npy * totals['exceedence_prob'])
            totals = totals.sortby('return_period')
            return totals

    tree = xr.DataTree()
    tree['ERA5'] = xr.DataTree(train_damages)
    tree['HazGAN']  = xr.DataTree(gener_damages)
    tree['Independent'] = xr.DataTree(indep_damages)
    tree = tree.map_over_datasets(calculate_total_return_periods, 'expected_damage')

    tree['Dependent'] = depen_damages
    tree['Dependent']['expected_damage'] = tree['Dependent']['expected_damage'].sum(dim=['lat', 'lon'])

    # clip return periods to a range of interest
    def truncate_rps(ds:xr.Dataset, minrp:float=1, maxrp:float=500) -> xr.Dataset:
        if len(ds.data_vars) > 0:
            ds = ds.where(ds['return_period'] >= minrp, drop=True)
            ds = ds.where(ds['return_period'] <= maxrp, drop=True)
        return ds

    tree = tree.map_over_datasets(truncate_rps)
    tree.to_netcdf(os.path.join(wd, "damage_scenarios.nc"))

# %%