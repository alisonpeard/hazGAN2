import numpy as np
import xarray as xr
from tqdm import tqdm
from typing import List, Union

from . import base
from .empirical import quantile
from .empirical import semiparametric_quantile


def invPIT(
        u:np.ndarray,
        x:np.ndarray,
        theta:np.ndarray=None,
        domain:Union[str, None]=None,
        distns:List[str]=["weibull", "genpareto", "genpareto"],
        two_tailed:List[bool]=[False, False, False]
    ) -> np.ndarray:
    """
    Transform uniform marginals to original distributions via
    inverse interpolation of empirical CDF.
    
    Parameters
    ----------
    u : np.ndarray
        Uniform marginals with shape [n, h, w, c] or [n, h * w, c]
    x : np.ndarray
        Original data for quantile calculation
    theta : np.ndarray, optional (default = None)
        Parameters of fitted Generalized Pareto Distribution (GPD)
    gumbel_margins : bool, optional (default = False)
        Whether to apply inverse Gumbel transform
    distribution : list, optional (default = ["weibull", "genpareto", "genpareto"])
        Distributions to use for quantile calculation. 

    Returns
    -------
    np.ndarray
        Transformed marginals with same shape as input u
    """
    u = u if domain is None else getattr(base, f"inv_{domain}")(u)
    
    original_shape = u.shape

    if u.ndim == 4:
        n, h, w, c = u.shape
        hw = h * w
        u = u.reshape(n, hw, c)
        x = x.reshape(len(x), hw, c)
        if theta is not None:
            theta = theta.reshape(hw, 6, c)
            theta = theta.transpose(1, 0, 2)
    elif u.ndim == 3:
        n, hw, c = u.shape
        if theta is not None:
            theta = theta.transpose(1, 0, 2)
    else:
        raise ValueError(
            "Uniform marginals must have dimensions [n, h, w, c] or [n, h * w, c]."
            )

    def transform(x, u, theta, i, c):
        """vectorised numpy transform."""
        x_i = x[:, i, c]
        u_i = u[:, i, c]
        theta_i = theta[:, i, c] if theta is not None else None
        distn = distns[c]
        tails = two_tailed[c]
        return (
            semiparametric_quantile(
                x_i, theta_i, distn, tails
                )(u_i)
            if theta is not None
            else quantile(x_i)(u_i)
        )
        
    quantiles = np.array([
        transform(x, u, theta, i, channel)
        for i in tqdm(range(hw)) for channel in range(c) 
    ])

    quantiles = quantiles.T
    quantiles = quantiles.reshape(*original_shape)

    return quantiles


def invPITDataset(
        ds:xr.Dataset,
        theta_var:str="params",
        u_var:str="uniform",
        x_var:str="anomaly",
        domain:Union[str, None]=None,
        distns:List[str]=["weibull", "genpareto", "genpareto"],
        two_tailed:List[bool]=[False, False, False]
        ) -> xr.DataArray:
    """
    Wrapper of invPIT for xarray.Dataset.
    """
    u = ds[u_var].values
    x = ds[x_var].values
    theta = ds[theta_var].values if theta_var in ds else None

    x_inv = invPIT(u, x, theta, domain, distns=distns, two_tailed=two_tailed)
    x_inv = xr.DataArray(x_inv, dims=ds[x_var].dims, coords=ds[x_var].coords)

    return x_inv


