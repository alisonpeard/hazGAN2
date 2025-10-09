# %%
import numpy as np
from .base import *

def chi_rmse(real, fake):
    c = real.shape[1]
    chi_diff = np.zeros(c)
    for channel in range(c):
        real_channel = real[:, channel, :, :]
        fake_channel = fake[:, channel, :, :]
        chi_diff[channel] = chi_diff_1d(real_channel, fake_channel)
    return np.mean(chi_diff)

def chi_diff_1d(real, fake):
    ecs_real = pairwise_extremal_coeffs(real)
    ecs_fake = pairwise_extremal_coeffs(fake)
    diff = ecs_real - ecs_fake
    return np.sqrt(np.mean(np.square(diff)))


def pairwise_extremal_coeffs(uniform):
    """Calculate extremal coefficients for each pair of pixels across single channel."""
    assert (
        len(uniform.shape) == 3
            ), "Function all_extremal_coeffs fors for 3-dimensional tensors only."
    n, h, w = uniform.shape
    uniform = np.reshape(uniform, (n, h * w))
    frechet = inverted_frechet(uniform)
    print("frechet shape:", frechet.shape)
    minima = minner_product(frechet.T, frechet)
    n = float(n)
    minima = minima.astype(float)
    ecs = np.zeros_like(minima)
    ecs = np.divide(n, minima, out=ecs, where=minima != 0, dtype=float)
    return ecs


def minner_product(a, b, batch_size=100):
    """Use broadcasting with batching to get sum of pairwise minima."""
    # a is (4096, 1248), b is (1248, 4096)
    n = a.shape[0]
    assert n == b.shape[-1], "First and last dimensions must match."
    result = np.zeros((n, n))
    
    for i in range(0, n, batch_size):
        batch_end = min(i + batch_size, n)
        batch_a = a[i:batch_end]  # Shape: (batch_size, 1248)
        
        batch_result = np.sum(
            np.minimum(
                np.expand_dims(batch_a, axis=-1),      # Shape: (batch_size, 1248, 1)
                np.expand_dims(b, axis=0)              # Shape: (1, 1248, 4096)
            ),
            axis=1
        )
        result[i:batch_end, :] = batch_result
    
    return result


def maxer_product(a, b):
    "Use broadcasting to get sum of pairwise maxima."
    x = np.sum(
            np.maximum(
                np.expand_dims(a, axis=-1),
                np.expand_dims(b, axis=0)),
            axis=1
        )
    return x


def test_minner_product():
    x = np.array([[1, 2], [1, 1]])
    assert np.array_equal(minner_product(x.T, x), np.array([[2, 2], [2, 3]])), f"{minner_product(x.T, x)}"


def get_extremal_coeffs_nd(marginals, sample_inds):
    """Calculate extremal coefficients across D-dimensional uniform data."""
    n, h, w, d = marginals.shape
    data = marginals.reshape(n, h * w, d)
    data = data[:, sample_inds, :]
    frechet = inverted_frechet(data)
    ecs = {}
    for i in range(len(sample_inds)):
        ecs[sample_inds[i]] = raw_extremal_coeff_nd(frechet[:, i, :])
    return ecs


def raw_extremal_coeff_nd(frechets):
    n, d = frechets.shape
    minima = np.min(frechets, axis=1)  # minimum for each row
    minima = np.sum(minima)
    if minima > 0:
        theta = n / minima
    else:
        print("Warning: all zeros in minima array.")
        theta = d
    return float(theta)


def _tail_dependence_coeff(u, v):
    """
    Classical tail dependence coefficient λ for upper tail dependence.
    Uses the Reiss & Thomas (2007) estimator.

    Args:
        u, v: 1D arrays of uniform marginals
    Returns:
        λ: tail dependence coefficient

    Refs:
        Reiss, R.-D. and Thomas, M. (2007) Statistical Analysis of Extreme Values,
        Birkhäuser, 3rd edition, Eq (2.62)
        https://rdrr.io/cran/extRemes/man/taildep.html
    """
    n = len(u)
    thresholds = np.arange(0.8, 0.99, 0.01)  # Multiple thresholds
    lambdas = []
    
    for t in thresholds:
        both_exceed = (u > t) & (v > t)
        
        # Reiss & Thomas (2007) formula: chi_hat = joint_exceedances / (n * (1-t))
        chi_t = np.sum(both_exceed) / (n * (1 - t))
        lambdas.append(chi_t)

    return np.mean(lambdas) if lambdas else 0