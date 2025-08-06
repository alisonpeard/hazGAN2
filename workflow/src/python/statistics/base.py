"""
Functions for transforming data to other distributions using empirical cdf.
"""
import numpy as np
from warnings import warn


def frechet(uniform):
    """Fréchet"""
    return -1 / np.log(uniform)


def inverted_frechet(uniform):
    """Inverted Fréchet RV is also exponentially distributed."""
    return -np.log(uniform)


def exp(uniform):
    """Exponential"""
    return -np.log(1 - uniform)


def inv_exp(uniform):
    """Inverse exponential"""
    return 1 - np.exp(-uniform)


def gumbel(uniform):
    """uniform -> Gumbel(0, 1)"""
    maxval = np.max(uniform)
    if maxval == 1:
        warn("Values == 1 found, scaling by 1e-6")
        uniform *= 1 - 1e-6
    if maxval > 1:
        raise ValueError(f"Some uniform > 1 ({maxval})")
    return -np.log(-np.log(uniform))


def inv_gumbel(x):
    """Gumbel(0, 1) -> uniform"""
    return np.exp(-np.exp(-x))


def laplace(uniform, mu=0, b=1):
    """uniform -> Laplace(mu, b) (quantile function)"""
    maxval = np.max(uniform)
    if maxval == 1:
        warn("Values == 1 found, scaling by 1e-6")
        uniform *= 1 - 1e-6
    if maxval > 1:
        raise ValueError(f"Some uniform > 1 ({maxval})")
    
    return np.where(
        uniform <= 0.5, 
        mu + b * np.log(2 * uniform),
        mu - b * np.log(2 - 2 * uniform)
        )


def inv_laplace(x, mu=0, b=1):
    """Laplace(mu, b) -> uniform (CDF function)."""
    return np.where(
        x <= mu,
        0.5 * np.exp((x - mu) / b),
        1 - 0.5 * np.exp(-(x - mu) / b)
    )
