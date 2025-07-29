"""
TODO: This code is super bloated, it can be made way more concise.

No need to turn into a data frame?
"""
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import logging

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.stats


TAIL = "upper" # temporarily hardcoded, should be a parameter


def main(input, output, params):
    EVENTS  = input.events
    FIGa = output.figa
    FIGb = output.figb
    FIGc = output.figc
    FIELDS  = params.fields
    PCRIT   = params.pcrit
    CMAP    = params.cmap
    FIGURES = [FIGa, FIGb, FIGc]


    def pk(var, tail=TAIL):
        return f"pk_{tail}_{var}"

    def thresh(var, tail=TAIL):
        return f"thresh_{tail}_{var}"

    def scale(var, tail=TAIL):
        return f"scale_{tail}_{var}"

    def shape(var, tail=TAIL):
        return f"shape_{tail}_{var}"


    # load GPD-fitted data   
    df = pd.read_parquet(EVENTS)
    df["lat"] = df["lat"].astype(float)
    df["lon"] = df["lon"].astype(float)
    df.columns = [col.replace(".", "_") for col in df.columns]
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(
        df["lon"],
        df["lat"],
        crs="EPSG:4326"
    )).set_crs("EPSG:4326")

    # log statistics (max, min) for each (variable, parameter)
    for var in FIELDS.keys():
        logging.info(f"Statistics for {var}_{TAIL}:")
        logging.info(f"{pk(var)}: {gdf[f'{pk(var)}'].min()} - {gdf[f'{pk(var)}'].max()}")
        logging.info(f"{thresh(var)}: {gdf[f'{thresh(var)}'].min()} - {gdf[f'{thresh(var)}'].max()}")
        logging.info(f"{scale(var)}: {gdf[f'{scale(var)}'].min()} - {gdf[f'{scale(var)}'].max()}")
        logging.info(f"{shape(var)}: {gdf[f'{shape(var)}'].min()} - {gdf[f'{shape(var)}'].max()}")

    # set values to be floats
    for var in FIELDS.keys():
        logging.info(f"Filling nans for {var}_{TAIL} with 0")
        gdf[pk(var)] = gdf[pk(var)].fillna(0)
        gdf[thresh(var)] = gdf[thresh(var)].fillna(0)
        gdf[scale(var)] = gdf[scale(var)].fillna(0)
        gdf[shape(var)] = gdf[shape(var)].fillna(0)

        logging.info(f"Converting {var}_{TAIL} to float")
        gdf[pk(var)] = gdf[pk(var)].astype(float)
        gdf[thresh(var)] = gdf[thresh(var)].astype(float)
        gdf[scale(var)] = gdf[scale(var)].astype(float)
        gdf[shape(var)] = gdf[shape(var)].astype(float)

    # turn it into an xarray dataset
    logging.info("Converting GeoDataFrame to xarray Dataset")
    ds = gdf.set_index(['lat', 'lon', 'event']).to_xarray().isel(event=0)

    logging.info(f"Latitude dimensions: {gdf['lat'].nunique()=}")
    logging.info(f"Longitude dimensions: {gdf['lon'].nunique()=}")

    n = 0
    for var, info in FIELDS.items():
        logging.info(f"Plotting {var}_{TAIL}")
        fig, axs = plt.subplots(1, 4, figsize=(14, 3), sharex=True, sharey=True,
        subplot_kw={'projection': ccrs.PlateCarree()}
        )

        # [left, bottom, width, height]
        ax4 = fig.add_axes([0.675, 0.25, 0.15, 0.55])
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
        logging.info(f"Using distribution: {scipy_distn}")
        logging.info(f"Loaded distribution: {dist}")

        try:
            # make cbars horizontal
            cbar_kwargs = {"label": None, "shrink": 0.9, "orientation": "horizontal", "pad": 0.05}
            logging.info(f"Plotting {var}_{TAIL}")
            logging.info(f"Plotting {pk(var)}")
            ds[pk(var)].plot(ax=axs[0], cmap=p_cmap, vmin=PCRIT, cbar_kwargs=cbar_kwargs, vmax=1)
            logging.info(f"Plotting {thresh(var)}")
            ds[thresh(var)].plot(ax=axs[1], cmap=cmap, cbar_kwargs=cbar_kwargs)
            logging.info(f"Plotting {scale(var)}")
            ds[scale(var)].plot(ax=axs[2], cmap=cmap, cbar_kwargs=cbar_kwargs)
            logging.info(f"Plotting {shape(var)}")
            ds[shape(var)].plot(ax=axs[3], cmap=cmap, cbar_kwargs=cbar_kwargs) #, vmin=-0.81, vmax=0.28)
            logging.info("Plotting complete")

            logging.info("Adding geo-features")
            for ax in[axs[0], axs[1], axs[2], axs[3]]:
                ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color="k")
                gl = ax.gridlines(draw_labels=False, linewidth=0.5, color='white', alpha=0.1)
                gl.top_labels = False
                gl.right_labels = False
                gl.left_labels = False
                gl.bottom_labels = False

            # plot some densities
            logging.info(f"Plotting densities for {var}_{TAIL}")
            shapes_all = gdf[shape(var)].values
            percentiles = np.linspace(0.01, 0.99, 10)
            shapes =  gdf[f'{shape(var)}'].quantile(percentiles)

            vmin = min(shapes_all)
            vmax = max(shapes_all)
            norm = plt.Normalize(vmin, vmax)
            colors = [plt.get_cmap(cmap)(norm(value)) for value in shapes]
        
            # plot the density for three sample grid points
            for i, shape_i in enumerate(shapes):
                u = np.linspace(0.95, 0.999, 100)
                x = dist.ppf(u, shape_i)
                y = dist.pdf(x, shape_i)
                ax4.plot(x, y, label=f"ξ={shape_i:.2f}", color=colors[i])

                # axis cleanup
                def percentage_formatter(x, pos):
                    return f'{100 * x:.0f}%'
                
                ax4.set_xlabel("")
                ax4.set_ylabel("")
                ax4.yaxis.set_major_formatter(percentage_formatter)
                ax4.tick_params(direction='in')
                ax4.yaxis.set_label_position("right")
            
            plt.close(fig)
        
        except Exception as e:
            logging.error(f"Error plotting {var}_{TAIL}: {e}")
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
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

        # plt.tight_layout()
        fig.savefig(FIGURES[n], bbox_inches="tight", dpi=300)
        logging.info(f"Saved figure {FIGURES[n]} for {var}_{TAIL}")
        n += 1


if __name__ == "__main__":
    # set up logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    main(input, output, params)