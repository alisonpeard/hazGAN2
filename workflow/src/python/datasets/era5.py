import os

# add new relevant variables here 
long_names = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "msl": "mean_sea_level_pressure",
    "tp": "total_precipitation",
    "ssrd": "surface_solar_radiation_downwards", 
    "t2m": "2m_temperature",
    "utci": "universal_thermal_climate_index",
    "lai": "leaf_area_index_vegetation",
    "i10fg": "instantaneous_10m_wind_gust",
    "swvl1": "volumetric_soil_water_layer_1",
    "sro": "surface_runoff" # accumulated
}

def get_input_file_pattern(indir, field):
    """Parse the input file pattern for the given field."""
    return os.path.join(indir, long_names[field], 'nc', '*.nc')


def filter_files(files, year, antecedent_buffer_days=None):
    """Filter the files for the given year."""
    files_year = [f for f in files if str(year) in f]

    if antecedent_buffer_days:
        nyears = 1 + (antecedent_buffer_days // 365)
        for dyear in range(1, nyears + 1):
            antecedent_files = [f for f in files if str(eval(year) - dyear) in f]
        files_year.extend(antecedent_files)

    return files_year


def clip_to_bbox(ds, xmin, xmax, ymin, ymax):
    """Clip the dataset to the given bounding box."""
    return ds.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))


def preprocess(ds):
    """Extra handline to preprocess the dataset."""
    for coord in ["expver", "number"]:
        if coord in ds.dims:
            ds = ds.sel({coord: ds[coord][0]})
        if coord in ds.coords:
            ds = ds.reset_coords(coord, drop=True)
    return ds