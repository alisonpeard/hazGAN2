"""Functions for plotting spatial correlation structures.
"""
import numpy as np
import matplotlib.pyplot as plt

from .base import CMAP
from .base import makegrid
from ..statistics import pairwise_extremal_coeffs


def smith1990(array):
    array = array.astype(np.float16)
    return pairwise_extremal_coeffs(array)


def pearson(array):
    n, h, w = array.shape
    array = array.reshape(n, h * w)
    corrs = np.corrcoef(array.T)
    print(corrs.shape)
    return corrs


def plot(fake, train, func, field=0, figsize=1., cmap=CMAP, vmin=None, vmax=None,
         title="Untitled", cbar_label="", alpha:float=1e-4, **func_kws):
    train = train[..., field].copy()
    fake  = fake[..., field].copy()
    
    train_res = func(train)
    fake_res  = func(fake)

    if alpha is not None:
        vmin = vmin or np.nanquantile(np.concatenate([fake_res, train_res]), alpha)
        vmax = vmax or np.nanquantile(np.concatenate([fake_res, train_res]), 1-alpha)
    else:
        vmin = vmin or train_res.min()
        vmax = vmax or train_res.max()
    
    cmap = getattr(plt.cm, cmap)
    cmap.set_under(cmap(0))
    cmap.set_over(cmap(.99))

    fig, axs, cax = makegrid(1, 2, figsize=figsize, projection=None)

    im = axs[0].imshow(train_res, cmap=cmap)
    _  = axs[-1].imshow(fake_res, cmap=cmap)

    for ax in axs:
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)

    axs[0].set_title("ERA5", y=-0.15)
    axs[-1].set_title("HazGAN", y=-0.15)
    fig.suptitle(title, y=1.05)

    fig.colorbar(im, cax=cax, label=cbar_label)

    return fig