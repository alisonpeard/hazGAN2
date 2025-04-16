"""Simple script to construct a global ERA5 grid for making parameters.

Not included in snakemake workflow as doesn't need to be re-run. Grid file
(~2 MB) is included in the resources folder.

rescaling factors from https://doi.org/10.1016/j.rser.2019.04.065
turbine 3PLE parameters from Table 4 https://doi.org/10.1016/j.epsr.2017.10.028
"""
# %%
import xarray as xr
from datetime import datetime

PROJECT = "renewablesuk"
input_file = "/Users/alison/Downloads/db284e6db9f28be19af6109191395be2.nc"
output_file = f"../resources/params/{PROJECT}.nc" # poject name

solar_params = {'g_stc': 1000, 't_ref': 25, 'c_0': 4.3, 'c_1': 0.943, 'c_2': 0.028, 'c_3': -1.528,}
onshore_wind_params = {'h': 80, 'gamma': 0.143, 'alpha': 2014.2, 'beta': 0.6846, 'ws0': 9.0431}
offshore_wind_params = {'h': 120, 'gamma': 0.111, 'alpha': 7038.5, 'beta': 0.7615, 'ws0': 8.7436}
onshore_params = solar_params | onshore_wind_params
offshore_params = solar_params | offshore_wind_params

def assign_params(land, onshore=onshore_params, offshore=offshore_params):
    """Assign turbine parameters based on land mask"""
    params = xr.where(land > 0.5, onshore, offshore)
    return params

if __name__ == "__main__":
    ds = xr.open_dataset(input_file)
    ds = ds[['latitude', 'longitude', 'lsm']]
    ds = ds.isel(valid_time=0)
    ds = ds.drop(['number', 'expver', 'valid_time'])
    ds = ds.rename({'lsm': 'land'})

    params = assign_params(ds.land)

    time = datetime.now()
    time = time.strftime("%Y-%m-%d %H:%M:%S")
    params.attrs = ds.attrs
    params.attrs['title'] = f"ERA5 grid for {PROJECT} project"
    params.attrs['created'] = time
    params.attrs['author'] = 'alison.peard@ouce.ox.ac.uk'

    # save a compressed netCDF file
    params = params.astype('float32')
    params.to_netcdf(output_file, encoding={
        var_name: {'zlib': True, 'complevel': 4} 
        for var_name in params.data_vars
    })

