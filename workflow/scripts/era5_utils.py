import numpy as np

long_names = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "msl": "mean_sea_level_pressure",
    "tp": "total_precipitation",
    "t2m": "2m_temperature",
    "ssrd": "surface_solar_radiation_downwards"
}

# variable creation functions
# data[field] = getattr(era5, func)(*[data[i] for i in infields])
def identity(x):
    return x

def wind_speed(u, v):
    """
    Compute the l2-norm of two arrays.
    """
    return np.sqrt(u**2 + v**2)



