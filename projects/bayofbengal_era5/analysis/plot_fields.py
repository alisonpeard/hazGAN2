# %% Damage probability field samples
import os
import xarray as xr
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# turn on spines
plt.rcParams['axes.spines.top'] = True
plt.rcParams['axes.spines.right'] = True


def damagefield(
        tree:xr.DataTree, label:str, ds:xr.Dataset, ds_var:str, 
        axs:plt.Axes, rp:int=10, add_colorbar:bool=False
        ) -> None:
    ref = tree[label].to_dataset()
    idx = (abs(ref['return_period'] - rp)).argmin().item()
    rp = ref['return_period'].isel(sample=idx).item()
    im1 = ds.isel(sample=idx, field=0)[ds_var].plot.contourf(
        ax=axs[0], cmap="viridis", levels=10, add_colorbar=False, vmin=0, vmax=32
        )
    im2 = ds.isel(sample=idx, field=1)[ds_var].plot.contourf(
        ax=axs[1], cmap="PuBu", levels=10, add_colorbar=False, vmin=0, vmax=0.66
        )
    im3 = ds.isel(sample=idx)['damage_prob'].plot.contourf(
        ax=axs[2], cmap="YlOrRd", levels=10, add_colorbar=False, vmin=0, vmax=0.9
        )
    axs[0].set_title(f"{label}\n(1-in-{rp:.0f})")
    axs[1].set_title("")
    axs[2].set_title("")

    if add_colorbar:
        # store the image references for later colorbar creation
        return im1, im2, im3
    
    


if __name__ == "__main__":
    
    wd = os.path.join("..", "results", "mangroves")
    tree = xr.open_dataset(os.path.join(wd, "damage_scenarios.nc"))
    train_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "train.nc"))
    gener_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "gener.nc"))
    depen_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "depen.nc"))
    indep_damages = xr.open_dataset(os.path.join(wd, "damage_fields", "indep.nc"))

    # Create the figure with the appropriate layout
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(3, 5, width_ratios=[1, 1, 1, 1, 0.05], wspace=0.1, hspace=0.1)

    # Create axes with projections
    axs = []
    for row in range(3):
        row_axes = []
        for col in range(4):
            ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
            row_axes.append(ax)
        axs.append(row_axes)
    axs = plt.np.array(axs)

    # Create separate axes for colorbars
    cbar_ax1 = fig.add_subplot(gs[0, 4])
    cbar_ax2 = fig.add_subplot(gs[1, 4])
    cbar_ax3 = fig.add_subplot(gs[2, 4])

    # Call damagefield for the first three columns without colorbars
    damagefield(tree, 'ERA5', train_damages, 'train', axs[:, 0])
    damagefield(tree, 'Dependent', depen_damages, 'dependent', axs[:, 1])
    damagefield(tree, 'Independent', indep_damages, 'independent', axs[:, 2])

    # Call the last column with colorbar flag set to True to get the image references
    im1, im2, im3 = damagefield(tree, 'HazGAN', gener_damages, 'fake', axs[:, 3], add_colorbar=True)

    # Add colorbars using the image references from the last column
    plt.colorbar(im1, cax=cbar_ax1, orientation='vertical', label='')
    plt.colorbar(im2, cax=cbar_ax2, orientation='vertical', label='')
    plt.colorbar(im3, cax=cbar_ax3, orientation='vertical', label='')

    # format colorbars to have percentage labels
    cbar_ax1.yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter("{x:.2f}"))
    cbar_ax2.yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter("{x:.2f}"))
    cbar_ax3.yaxis.set_major_formatter(plt.matplotlib.ticker.PercentFormatter(1, 0))

    # Set labels and features for all axes
    for ax in axs.flat:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.2)

    axs[0, 0].set_ylabel(r"Wind speed (ms$^{-1}$)")
    axs[1, 0].set_ylabel("Precipitation (m)")
    axs[2, 0].set_ylabel("Damage probability")

    # Add titles to the top row
    axs[-1, 0].set_title("ERA5\n(1-in-10)", y=-.4)
    axs[-1, 1].set_title("Dependent\n(1-in-10)", y=-.4)
    axs[-1, 2].set_title("Independent\n(1-in-10)", y=-.4)
    axs[-1, 3].set_title("HazGAN\n(1-in-10)", y=-.4)

    plt.tight_layout()
    plt.show()

    os.makedirs('../results/figures/mangroves', exist_ok=True)
    fig.savefig("../results/figures/mangroves/damage_fields.png", dpi=300, bbox_inches='tight', transparent=True)

# %%