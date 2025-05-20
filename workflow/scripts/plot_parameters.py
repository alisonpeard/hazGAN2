import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.stats
import logging

if __name__ == "__main__":
    DATA_ALL= snakemake.input.data_all
    EVENTS  = snakemake.input.events
    FIGa = snakemake.output.figa
    FIGb = snakemake.output.figb
    FIGc = snakemake.output.figc
    FIELDS  = snakemake.params.fields
    PCRIT   = snakemake.params.pcrit
    CMAP    = snakemake.params.cmap

    FIGURES = [FIGa, FIGb, FIGc]

    # set up logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    # load coordinates
    coords = xr.open_dataset(DATA_ALL)
    coords = coords[["lat","lon"]].to_dataframe().reset_index()
    coords = gpd.GeoDataFrame(
        coords, geometry=gpd.points_from_xy(coords['lon'], coords['lat'])
        ).set_crs("EPSG:4326")

    # load GPD-fitted data   
    df = pd.read_parquet(EVENTS)
    df["lat"] = df["lat"].astype(float).round(6) # about 11 cm
    df["lon"] = df["lon"].astype(float).round(6)
    coords["lat"] = coords["lat"].astype(float).round(6)
    coords["lon"] = coords["lon"].astype(float).round(6)
    df = df.merge(coords, on=["lat", "lon"])
    df.columns = [col.replace(".", "_") for col in df.columns]
    gdf = gpd.GeoDataFrame(df, geometry="geometry").set_crs("EPSG:4326")

    # log statistics (max, min) for each (variable, parameter)
    for var in FIELDS.keys():
        logging.info(f"Statistics for {var}:")
        logging.info(f"pk_{var}: {gdf[f'pk_{var}'].min()} - {gdf[f'pk_{var}'].max()}")
        logging.info(f"thresh_{var}: {gdf[f'thresh_{var}'].min()} - {gdf[f'thresh_{var}'].max()}")
        logging.info(f"scale_{var}: {gdf[f'scale_{var}'].min()} - {gdf[f'scale_{var}'].max()}")
        logging.info(f"shape_{var}: {gdf[f'shape_{var}'].min()} - {gdf[f'shape_{var}'].max()}")

    # set values to be floats
    for var in FIELDS.keys():
        logging.info(f"Filling nans for {var} with 0")
        gdf[f"pk_{var}"] = gdf[f"pk_{var}"].fillna(0)
        gdf[f"thresh_{var}"] = gdf[f"thresh_{var}"].fillna(0)
        gdf[f"scale_{var}"] = gdf[f"scale_{var}"].fillna(0)
        gdf[f"shape_{var}"] = gdf[f"shape_{var}"].fillna(0)

        logging.info(f"Converting {var} to float")
        gdf[f"pk_{var}"] = gdf[f"pk_{var}"].astype(float)
        gdf[f"thresh_{var}"] = gdf[f"thresh_{var}"].astype(float)
        gdf[f"scale_{var}"] = gdf[f"scale_{var}"].astype(float)
        gdf[f"shape_{var}"] = gdf[f"shape_{var}"].astype(float)

    # turn it into an xarray dataset
    ds = gdf.set_index(['lat', 'lon', 'event']).to_xarray().isel(event=0)

    n = 0
    for var, info in FIELDS.items():
        fig, axs = plt.subplots(1, 4, figsize=(14, 3), sharex=True, sharey=True,
        subplot_kw={'projection': ccrs.PlateCarree()}
        )

        # [left, bottom, width, height]
        # ax4 = fig.add_axes([0.825, 0.1, 0.1, 0.8])
        ax4 = fig.add_axes([0.675, 0.25, 0.15, 0.55])
        # ax4 = fig.add_axes([0.8, 0, 0.2, 1])
        plt.subplots_adjust(right=0.65) 

        cmap = CMAP
        p_cmap = plt.get_cmap(cmap)
        p_cmap.set_under("crimson")

        # plot the density for three sample grid points
        if info["distn"] == "weibull":
            scipy_distn = "weibull_min"
        else:
            scipy_distn = info["distn"]
        dist = getattr(scipy.stats, scipy_distn)

        try:
            # make cbars horizontal
            cbar_kwargs = {"label": None, "shrink": 0.7, "orientation": "horizontal", "pad": 0.05}
            ds[f"pk_{var}"].plot(ax=axs[0], cmap=p_cmap, vmin=PCRIT, cbar_kwargs=cbar_kwargs)
            ds[f"thresh_{var}"].plot(ax=axs[1], cmap=cmap, cbar_kwargs=cbar_kwargs)
            ds[f"scale_{var}"].plot(ax=axs[2], cmap=cmap, cbar_kwargs=cbar_kwargs)
            ds[f"shape_{var}"].plot(ax=axs[3], cmap=cmap, cbar_kwargs=cbar_kwargs) #, vmin=-0.81, vmax=0.28)

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
        
            # plot the density for three sample grid points
            for i, shape in enumerate(shapes):
                u = np.linspace(0.95, 0.999, 100)
                x = dist.ppf(u, shape) #, loc=loc, scale=scale)
                y = dist.pdf(x, shape) #, loc=loc, scale=scale)
                ax4.plot(x, y, label=f"ξ={shape:.2f}", color=colors[i])

                # axis cleanup
                def percentage_formatter(x, pos):
                    return f'{100 * x:.0f}%'
                
                ax4.set_xlabel("")
                ax4.set_ylabel("")
                ax4.yaxis.set_major_formatter(percentage_formatter)
                ax4.tick_params(direction='in')
                ax4.yaxis.set_label_position("right")
        
        except Exception as e:
            logging.error(f"Error plotting {var}: {e}")
            continue

        # # add colorbar for shape
        # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        # sm.set_array([])
        # cbar = fig.colorbar(sm, ax=ax4, orientation="horizontal", pad=0.05)
        # cbar.set_label(None)
        # cbar.set_label("ξ")
        # cbar.ax.yaxis.set_label_position('left')
        # cbar.ax.set_ylabel('ξ', rotation=0, labelpad=15)

        # turn off top and right spines for ax4
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)

        # ax4.legend()
        axs[0].set_title("H0: X~{}".format(info["distn"].capitalize()))
        axs[1].set_title("μ")
        axs[2].set_title("σ")
        axs[3].set_title("ξ")

        for ax in axs[:-1].ravel():
            # ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

        # plt.tight_layout()
        fig.savefig(FIGURES[n], bbox_inches="tight", dpi=300)
        n += 1