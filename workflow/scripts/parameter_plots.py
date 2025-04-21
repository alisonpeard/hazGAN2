import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.stats

if __name__ == "__main__":
    DATA_ALL= snakemake.input.data_all
    EVENTS  = snakemake.input.events
    FIGURES = snakemake.output.figures
    FIELDS  = snakemake.params.fields
    PCRIT   = snakemake.params.pcrit
    CMAP    = snakemake.params.cmap

    # ds = xr.open_dataset(DATA)

    # load coordinates
    coords = xr.open_dataset(DATA_ALL)
    coords = coords['grid'].to_dataframe().reset_index()
    coords = gpd.GeoDataFrame(
        coords, geometry=gpd.points_from_xy(coords['lon'], coords['lat'])
        ).set_crs("EPSG:4326")

    # load GPD-fitted data   
    df = pd.read_parquet(EVENTS)
    df = df.merge(coords, on="grid")
    df.columns = [col.replace(".", "_") for col in df.columns]
    gdf = gpd.GeoDataFrame(df, geometry="geometry").set_crs("EPSG:4326")

    # turn it into an xarray dataset
    ds = gdf.set_index(['lat', 'lon', 'event']).to_xarray().isel(event=0)

    n = 0
    for var, info in FIELDS.items():
        # check variables aren't all NA
        if ds[f"shape_{var}"].isnull().all():
            print(f"Skipping {var} because all values are NaN")
            n += 1
            continue

        fig, axs = plt.subplots(1, 4, figsize=(16, 3), sharex=True, sharey=True,
        subplot_kw={'projection': ccrs.PlateCarree()})

        ax4 = fig.add_axes([0.825, 0.1, 0.15, 0.8])   # [left, bottom, width, height]
        plt.tight_layout()
        plt.subplots_adjust(right=0.8) 

        cmap = CMAP
        p_cmap = plt.get_cmap(cmap)
        p_cmap.set_under("crimson")

        ds[f"p_{var}"].plot(ax=axs[0], cmap=p_cmap, vmin=PCRIT, cbar_kwargs={'label': None})
        ds[f"thresh_{var}"].plot(ax=axs[1], cmap=cmap, cbar_kwargs={'label': None})
        ds[f"scale_{var}"].plot(ax=axs[2], cmap=cmap, cbar_kwargs={'label': None})
        ds[f"shape_{var}"].plot(ax=axs[3], cmap=cmap, add_colorbar=False) #, vmin=-0.81, vmax=0.28)

        # plot the density for three sample grid points
        if info["distn"] == "weibull":
            scipy_distn = "weibull_min"
        else:
            scipy_distn = info["distn"]
        dist = getattr(scipy.stats, scipy_distn)
        # if var == 'u10':
        #     dist = getattr(scipy.stats, 'weibull_min')
        # else:
        #     dist = getattr(scipy.stats, 'genpareto')

        # plot some densities
        shapes_all = gdf[f'shape_{var}'].values
        percentiles = np.linspace(0.01, 0.99, 10)
        shapes =  gdf[f'shape_{var}'].quantile(percentiles)
        loc    = gdf[f'thresh_{var}'].mean()
        scale  = gdf[f'scale_{var}'].mean()

        vmin = min(shapes_all)
        vmax = max(shapes_all)
        norm = plt.Normalize(vmin, vmax)
        colors = [plt.get_cmap(cmap)(norm(value)) for value in shapes]
        
        for i, shape in enumerate(shapes):
            u = np.linspace(0.95, 0.999, 100)
            x = dist.ppf(u, shape) #, loc=loc, scale=scale)
            y = dist.pdf(x, shape) #, loc=loc, scale=scale)
            ax4.plot(x, y, label=f"ξ={shape:.2f}", color=colors[i])

            # axis cleanup
            def percentage_formatter(x, pos):
                return f'{100 * x:.0f}%'  # Multiply by 100 and add % sign
            
            ax4.set_xlabel("")
            ax4.set_ylabel("")
            ax4.yaxis.set_major_formatter(percentage_formatter)
            ax4.tick_params(direction='in')
            ax4.yaxis.set_label_position("right")
            # ax4.tick_params(axis='y', pad=-40)

        # add colorbar for shape
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax4)
        # set label to none
        cbar.set_label(None)
        cbar.set_label("ξ")
        # put cbar label on LEFT of it
        cbar.ax.yaxis.set_label_position('left')
        # rotate the cbar label to be UPRIGHT
        cbar.ax.set_ylabel('ξ', rotation=0, labelpad=15)
        # cbar.set_ticks(shapes)
        # cbar.set_ticklabels([f"{s:.2f}" for s in shapes])

        # ax4.legend()
        axs[0].set_title("H0: X~{}".format(info["distn"].capitalize()))
        # if var == 'u10':
        #     axs[0].set_title("H₀: X~Weibull(ξ,μ,σ)")
        # else:
        #     axs[0].set_title("H₀: X~GPD(ξ,μ,σ)")

        axs[1].set_title("μ")
        axs[2].set_title("σ")
        axs[3].set_title("ξ")

        for ax in axs[:-1].ravel():
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

        # fig.suptitle(f"Fit for ERA5 {var.upper()}, n = {gdf['event'].nunique()}")
        # print(gdf[gdf[f"pk_{var}"] < pcrit]["grid"].nunique(), "significant p-values")
        fig.savefig(FIGURES[n], bbox_inches='tight', dpi=300)
        n += 1