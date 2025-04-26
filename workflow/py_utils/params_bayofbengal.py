"""Simple script to construct a global ERA5 grid for making parameters.

Download ERA5 data from CCDS.

Not included in snakemake workflow as doesn't need to be re-run. Grid file
(~2 MB) is included in the resources folder.
"""
# %%
import xarray as xr

input_file = "/Users/alison/Downloads/db284e6db9f28be19af6109191395be2.nc"
output_file = "../resources/params/bayofbengal.nc"

if __name__ == "__main__":
    params = xr.open_dataset(input_file)
    params = params[['latitude', 'longitude', 'lsm']]
    params = params.isel(valid_time=0)
    params = params.drop(['number', 'expver', 'valid_time'])
    params = params.rename_vars({"lsm": "sample"})
    params.sample.attrs['long_name'] = "Sample parameter"
    params.sample.attrs['units'] = "no units"
    params.attrs['info'] = "Sample parameter file with ERA5 grid. Made 16-04-2025."
    
    # save a compressed netCDF file
    params = params.astype('float32')
    params.to_netcdf(output_file, encoding={
        var_name: {'zlib': True, 'complevel': 4} 
        for var_name in params.data_vars
    })

# %%