
import os
from PIL import Image
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import pandas as pd
import xarray as xr
from hazGAN.statistics import invPIT, gumbel


__all__ = ['plot', 'load_samples', 'yflip']


CMAP  = "Spectral_r"
EPS   = 1e-6


def yflip(array: np.ndarray, ax=1) -> np.ndarray:
    if array is not None:
        return np.flip(array, axis=ax)
        

def load_samples(samples_dir, data_dir, model, threshold=None, sampletype='samples'):
    """Load and process samples and training data for visualisation"""
    # load samples - - - - - - - - --- - - - -- - -- - - -- -- - - -- - - - -
    samples_path = os.path.join(samples_dir, model, "results", sampletype)
    samples_list = glob(os.path.join(samples_path, "seed*.png"))
    samples_list = sorted(samples_list)
    print(f"Found {len(samples_list)} samples in {samples_path}")

    samples = []
    for png in samples_list:
        img = Image.open(png)
        samples.append(np.array(img))
    samples = np.array(samples).astype(float)
    samples /= 255.
    print(f"Loaded {samples.shape} samples")

    # load gumbel scaling statistics
    stats_file = os.path.join(data_dir, "images", "gumbel", "image_stats.npz")
    stats = np.load(stats_file)
    image_minima = stats['min']
    image_maxima = stats['max']
    n            = stats['n']
    image_range  = image_maxima - image_minima

    print(f"Loaded image statistics with shape {image_minima.shape}")

    # rescale images 
    Warning("Using new rescaling from Wed 22 January 2025")
    samples = (samples * (n + 1) - 1) / (n - 1) * image_range + image_minima

    # order samples
    sample_maxima = samples[..., 0].max(axis=(1,2))
    sample_order  = np.argsort(sample_maxima)[::-1]
    samples        = samples[sample_order]

    # check sample distribution (Gumbel)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.hist(samples.ravel(), bins=50, color='lightgrey', edgecolor='k', density=True);
    ax.set_xlabel("Pixel value");
    ax.set_ylabel("Density");
    ax.set_title("Histogram of all Gumbel(0, 1) samples") # this should be uniform

    # load training data - - - - - - - - --- - - - -- - -- - - -- -- - - -- - - - -
    data   = xr.open_dataset(os.path.join(data_dir, "data.nc"))
    data['maxwind'] = data.sel(field='u10')['anomaly'].max(dim=['lat', 'lon'])
    trainmask = data.where(data['time.year'] != 2021, drop=True).time
    validmask = data.where(data['time.year'] == 2021, drop=True).time
    valid = data.sel(time=validmask)
    data  = data.sel(time=trainmask)
    if threshold is not None:
        tmask  = data.where(data['maxwind'] >= 15., drop=True).time
        data   = data.sel(time=tmask)
    data   = data.sortby('maxwind', ascending=False)
    x      = data.anomaly.values
    u      = data.uniform.values
    params = data.params.values
    print(f"Loaded {params.shape} parameters")

    if threshold is not None:
        tmask = valid.where(valid['maxwind'] >= 15., drop=True).time
        valid = valid.sel(time=tmask)
    valid = valid.sortby('maxwind', ascending=False)
    x_valid = valid.anomaly.values
    u_valid = valid.uniform.values

    # order training samples by x
    x_maxima       = x[..., 0].max(axis=(1, 2))
    x_order        = np.argsort(x_maxima)[::-1]
    u              = u[x_order]
    x              = x[x_order]

    if (u >= 1).sum() > 0:
        print("Some data.nc is greater than 1")
        invalid_umask = (u >= 1).astype(bool)
        u *= (1-EPS)
    else:
        invalid_umask = None

    if (u_valid >= 1.).sum() > 0:
        print("Some valid data.nc is greater than 1")
        invalid_valid_umask = (u_valid >= 1).astype(bool)
        u_valid *= (1-EPS)
    else:
        invalid_valid_umask = None

    # transformations
    x_gumbel        = gumbel(u)
    valid_gumbel    = gumbel(u_valid)
    samples_uniform = np.exp(-np.exp(-samples))

    # remove [0,1] values
    if (samples_uniform >=1 ).sum() > 0:
        Warning("Some uniform samples are greater than 1")
        invalid_mask = (samples_uniform>=1).astype(bool)
        samples *= (1-EPS)
    else:
        invalid_mask = None

    if np.nanmin(samples_uniform) <= 0.:
        Warning("Some uniform samples == 0")
        samples_uniform = np.clip(samples_uniform, EPS, float('inf'))

    # get samples into data space
    samples_temp = np.flip(samples_uniform, axis=1)
    samples_x = invPIT(samples_temp, x, params)
    samples_x = np.flip(samples_x, axis=1)
    del samples_temp

    # reorder samples in x space
    sample_maxima   = samples_x[..., 0].max(axis=(1,2))
    sample_order    = np.argsort(sample_maxima)[::-1]
    samples_x       = samples_x[sample_order]
    samples_uniform = samples_uniform[sample_order]
    samples         = samples[sample_order]

    # negate MSLP
    samples_x[..., 2] *= -1
    x[..., 2]         *= -1
    x_valid[..., 2]   *= -1

    # these will be upside-down in heatmaps and correctly orientated
    # in contour plots because y goes from 80 -> 95
    return {
        'samples': {
            'uniform': yflip(samples_uniform),
            'gumbel':  yflip(samples),
            'data':    yflip(samples_x),
            'mask':    yflip(invalid_mask)
        },
        'training': {
            'uniform': u, # yflip(u),
            'gumbel':  x_gumbel, # yflip(x_gumbel),
            'data':    x, # yflip(x),
            'mask':    invalid_mask #yflip(invalid_umask, 0)
        },
        'valid': {
            'uniform': u_valid, # yflip(u_valid),
            'gumbel':  valid_gumbel, # yflip(valid_gumbel),
            'data':    x_valid, # yflip(x_valid),
            'mask':    invalid_valid_umask # yflip(invalid_valid_umask, 0)
        }
    }


def plot(array, field, yflip=False, contours=False, mask=None, title='',
         exclude_mask=False, print_stats=True,
         standardise_colours=True,
         vmin=None, vmax=None,
         cmap=CMAP, levels=13):
    """

    """
    array = array.copy()
    if yflip:
        array = np.flip(array, axis=1)

    if exclude_mask and mask is not None:
        array[mask] = np.nan

    if standardise_colours:
        vmin = vmin or np.nanmin(array[..., field])
        vmax = vmax or np.nanmax(array[..., field])

    fig, axs = plt.subplots(8, 8, figsize=(16, 13),
                            sharex=True, sharey=True,
                            gridspec_kw={'hspace': 0., 'wspace': 0.})
    for i, ax in enumerate(axs.flat):
        if contours:
            levels = np.linspace(vmin, vmax, levels)
            im = ax.contourf(array[i, ..., field],
            cmap=cmap,
            levels=levels
            )
        else:
            im = ax.imshow(array[i, ..., field],
                           cmap=cmap,
                           vmin=vmin, vmax=vmax
                           )
        if mask is not None:
            ax.contour(mask[i, ..., field],
                       colors='k',
                       linewidths=1
                       )
        ax.label_outer()
    
    if standardise_colours:
        fig.colorbar(im, ax=list(axs.flat))
    fig.suptitle(title, y=.9)

    if print_stats:
        print(f"\nStatistics for {title}:\n----------------------")
        print(f"Min: {np.nanmin(array[..., field]):.4f}")
        print(f"Max: {np.nanmax(array[..., field]):.4f}")
        print(f"Mean: {np.nanmean(array[..., field]):.4f}")
        print(f"Std: {np.nanstd(array[..., field]):.4f}")