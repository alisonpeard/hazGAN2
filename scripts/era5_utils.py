import numpy as np

long_names = {
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
    "msl": "mean_sea_level_pressure",
    "tp": "total_precipitation",
}

def identity(a):
    return a

def l2norm(a, b):
    """
    Compute the l2-norm of two arrays.
    """
    return np.sqrt(a**2 + b**2)



