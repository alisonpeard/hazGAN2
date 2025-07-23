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
    """
    Parse the input pattern for the given field.
    """
    return os.path.join(indir, long_names[field], 'nc', '*.nc')


def filter_files(files, year):
    """
    Filter the files for the given year.
    """
    return [f for f in files if str(year) in f]


def clip_to_bbox(ds, xmin, xmax, ymin, ymax):
    """
    Clip the dataset to the given bounding box.
    """
    return ds.sel(longitude=slice(xmin, xmax), latitude=slice(ymax, ymin))