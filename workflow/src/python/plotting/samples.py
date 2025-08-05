"""Functions for plotting data samples.
"""
import numpy as np
from .base import CMAP
from .base import makegrid, contourmap, scalebar
from ..statistics import invPIT


def identity(array, *args, **kwargs):
    return array


def gumbel(array, axis=0):
    array = np.clip(array, 1e-6, 1 - 1e-6)
    return -np.log(-np.log(array))


def anomaly(array, reference, params):
    array = invPIT(array, reference, params)
    return array


def plot(fake, train, field=0, transform=None,
    vmin=None, vmax=None, cmap=CMAP, extent=None,
    cbar_label='', cbar_width=0.2, linewidth=.1, alpha=1e-4, alpha_vlim=True, 
    nrows=4, ncols=8, ndecimals=1, **transform_kws):
    """Plot training samples on top row and generated samples on bottom row."""

    assert extent is not None, "Extent must be provided for samples plots."

    transform = transform or identity
    fake  = transform(fake, **transform_kws)
    train = transform(train, **transform_kws)

    fake = fake[..., field].copy()
    train = train[..., field].copy()

    if alpha_vlim:
        vmin = vmin or np.nanquantile(np.concatenate([train]), alpha) # used to be [fake, train]
        vmax = vmax or np.nanquantile(np.concatenate([train]), 1-alpha) # used to be [fake, train]
    else:
        vmin = vmin or min(np.nanmin(fake), np.nanmin(train))
        vmax = vmax or max(np.nanmax(fake), np.nanmax(train))

    nrows = 4
    ncols = 8
    total = nrows * ncols
    midpoint = total // 2
    midrow  = nrows // 2

    fig, axs, cax = makegrid(nrows, ncols, cbar_width=cbar_width, figsize=1.)
    for i, ax in enumerate(axs.flat):
        if i < midpoint:
            contourmap(fake[i, ...], ax=ax,
                       extent=extent,
                       vmin=vmin, vmax=vmax,
                       cmap=cmap, linewidth=linewidth, ndecimals=ndecimals)
        if i >= midpoint:
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0 - 0.01, pos.width, pos.height])
            j = i - midpoint
            im = contourmap(train[j, ...], ax=ax,
                            extent=extent,
                            vmin=vmin, vmax=vmax,
                            cmap=cmap, linewidth=linewidth, ndecimals=ndecimals)

    axs[0, 0].set_ylabel("HazGAN", fontsize=18)
    axs[midrow, 0].set_ylabel("ERA5", fontsize=18)
    
    scalebar(axs[-1, -1])
    scalebar(axs[midrow-1, -1])

    fig.colorbar(im, cax=cax, label=cbar_label)
    return fig