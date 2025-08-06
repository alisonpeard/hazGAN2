"""Base plotting functions and settings.
"""
import numpy as np
import cmocean.cm as cmo
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.geodesic as geodesic
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar


plt.rcParams['font.size'] = 12
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['legend.frameon'] = False  # Frameless legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

FIGSIZE_ONECOL = (3.5, 2.16)
FIGSIZE_ONECOL_LARGE = (3.5, 3.5)
FIGSIZE_TWOCOL = (7.0, 4.33)
FIGSIZE_TWOCOL_WIDE = (7.0, 3.0)


CMAP   = "viridis"
LABELS = ["wind speed [m]", "total precipitation [m]", "mslp [Pa]"]
hist_kwargs = {"bins": 50, "color": "lightgrey", "edgecolor": "black", "alpha": 0.7}
save_kwargs = {"dpi": 300, "bbox_inches": "tight", "pad_inches": 0.1}


def eval_cmap_str(name:str):
    """Return a cmocean colormap."""
    if name.startswith("cmocean.cm."):
        try:
            name = name.split(".")[-1]
            return getattr(cmo, name)
        except Exception as e:
            print(f"Couldn't retrive colormap {name} from cmocean.cm. Reason: {e}. Defaulting to {CMAP}.")

    return plt.get_cmap(name)


def linspace(start, stop, num=50, ndecimals=1):
        """Return a linspace with up to ndecimals decimal places."""
        factor = 10 ** ndecimals
        levels = np.linspace(start * factor, stop * factor, num, dtype=int) / factor
        return sorted(set(levels))


def makegrid(rows, cols, figsize:float=1, fig=None,
             cbar_width=.05, cbar_pad=.1, projection=ccrs.PlateCarree()):    

    total_width = cols + cbar_width + cbar_pad
    height = rows
    scale   = 8 / max(total_width, height)
    figsize = (figsize * total_width * scale, figsize * height * scale)
    
    if fig is None:
        fig = plt.figure(figsize=figsize)
    else:
        tmp_fig = plt.figure(figsize=figsize)
        width, height = tmp_fig.get_size_inches()
        fig.set_size_inches(width, height)
        plt.close(tmp_fig)
    
    gs = fig.add_gridspec(rows, cols + 2, 
                         width_ratios=[1]*cols + [cbar_pad, cbar_width],
                         wspace=0,
                         hspace=0,
                         left=0.01, right=0.99,
                         top=0.99, bottom=0.01)
    
    axes = np.empty((rows, cols), dtype=object)
    for i in range(rows):
        for j in range(cols):
            axes[i,j] = fig.add_subplot(gs[i, j], projection=projection)
            axes[i,j].set_xticks([])
            axes[i,j].set_yticks([])
    
    cax = fig.add_subplot(gs[:, -1])
    
    return fig, axes.squeeze(), cax


def scalebar(ax, location='lower right', length_fraction=.6):
    # calculate dx
    xlim = ax.get_xlim()
    npixels = int(xlim[1] - xlim[0])
    extent = ax.get_extent() # [80, 95, 10, 25]

    points    = (extent[0], extent[2])
    endpoints = (extent[1], extent[2])

    geod = geodesic.Geodesic()
    dx = geod.inverse(points, endpoints)
    dx = dx[0][0] / npixels
    scalebar = ScaleBar(dx, loc=location, fixed_value=1000, fixed_units='km',
                        box_alpha=.8, frameon=True, box_color='white', color='k',
                        bbox_to_anchor=(.8, .05), bbox_transform=ax.transAxes,
                        width_fraction=.05)
    ax.add_artist(scalebar)


def heatmap(array, ax=None, extent=None, transform=ccrs.PlateCarree(),
            cmap=CMAP, vmin=None, vmax=None, title=False, linewidth=.5):
    """Plot a heatmap with the coastline."""
    assert extent is not None, "Extent must be provided for heatmap."
    h, w = array.shape
    array = np.flip(array, axis=0) # imshow is upside down
    ax = ax or plt.axes(projection=transform)
    im = ax.imshow(array, extent=extent, transform=transform, cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation='nearest')
    ax.add_feature(cfeature.COASTLINE, edgecolor='k', linewidth=linewidth)
    ax.set_extent(extent)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=16);
    return im


def contourmap(array, ax=None, extent=None, transform=ccrs.PlateCarree(),
            cmap=CMAP, vmin=None, vmax=None, levels=10, extend="both", title=False,
            linewidth=.5, ndecimals=1):
    """Plot a heatmap with the coastline."""
    assert extent is not None, "Extent must be provided for contourmap."
    h, w = array.shape
    ax = ax or plt.axes(projection=transform)

    vmin = vmin or array.min()
    vmax = vmax or array.max()

    levels = linspace(vmin, vmax, levels, ndecimals)
    im = ax.contourf(array, extent=extent, transform=transform, cmap=cmap, levels=levels,
                     extend=extend)
    ax.add_feature(cfeature.COASTLINE, edgecolor='k', linewidth=linewidth)
    ax.set_extent(extent)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=16);
    return im