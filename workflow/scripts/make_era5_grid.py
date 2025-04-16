"""Simple script to construct a global ERA5 grid for making parameters.

Not included in snakemake workflow as doesn't need to be re-run. Grid file
(~2 MB) is included in the resources folder.
"""
import xarray as xr

input_file = "/Users/alison/Downloads/9ec0bcd2d4351542cb403cb6f73f500d.nc"
output_file = "../resources/grids/era5.nc"

if __name__ == "__main__":
    params = xr.open_dataset(input_file)
    params = params[['latitude', 'longitude', 'u10']]
    params = params.isel(valid_time=0)
    params = params.drop(['number', 'expver', 'valid_time'])
    params = params.rename_vars({"u10": "sample"})
    params.sample.attrs['long_name'] = "Sample parameter"
    params.sample.attrs['units'] = "no units"
    params.attrs['info'] = "Sample parameter file with ERA5 grid. Made 16-04-2025."
    
    # save a compressed netCDF file
    params = params.astype('float32')
    params.to_netcdf(output_file, encoding={
        var_name: {'zlib': True, 'complevel': 4} 
        for var_name in params.data_vars
    })

