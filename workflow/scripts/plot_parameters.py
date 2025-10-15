"""
TODO: This code is super bloated, it can be made way more concise.

No need to turn into a data frame?
"""
import numpy as np
import pandas as pd
import geopandas as gpd
# import xarray as xr
from functools import partial
import logging

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.stats


def pk(var, tail):
    return f"pk_{tail}_{var}"

def thresh(var, tail):
    return f"thresh_{tail}_{var}"

def scale(var, tail):
    return f"scale_{tail}_{var}"

def shape(var, tail):
    return f"shape_{tail}_{var}"


def plot_parameters_for_data(
        ds, gdf, var, fields, tail, cmap, pcrit, outpath
        ) -> None:
    logging.info(f"Plotting {var}_{tail}")
    info = fields[var]

    pk_n     = partial(pk, tail=tail)
    thresh_n = partial(thresh, tail=tail)
    scale_n  = partial(scale, tail=tail)
    shape_n  = partial(shape, tail=tail)

    fig, axs = plt.subplots(1, 4, figsize=(14, 3), sharex=True, sharey=True,
    subplot_kw={'projection': ccrs.PlateCarree()}
    )

    # [left, bottom, width, height]
    ax4 = fig.add_axes([0.675, 0.25, 0.15, 0.55])
    plt.subplots_adjust(right=0.65) 

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
        logging.info(f"Plotting {var}_{tail}")
        logging.info(f"Plotting {pk_n(var)}")
        ds[pk_n(var)].plot(ax=axs[0], cmap=p_cmap, vmin=pcrit, cbar_kwargs=cbar_kwargs, vmax=1)
        logging.info(f"Plotting {thresh_n(var)}")
        ds[thresh_n(var)].plot(ax=axs[1], cmap=cmap, cbar_kwargs=cbar_kwargs)
        logging.info(f"Plotting {scale_n(var)}")
        ds[scale_n(var)].plot(ax=axs[2], cmap=cmap, cbar_kwargs=cbar_kwargs)
        logging.info(f"Plotting {shape_n(var)}")
        ds[shape_n(var)].plot(ax=axs[3], cmap=cmap, cbar_kwargs=cbar_kwargs) #, vmin=-0.81, vmax=0.28)
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
        logging.info(f"Plotting densities for {var}_{tail}")
        shapes_all = gdf[shape_n(var)].values
        percentiles = np.linspace(0.01, 0.99, 10)
        shapes =  gdf[f'{shape_n(var)}'].quantile(percentiles)

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
        logging.error(f"Error plotting {var}_{tail}: {e}")

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
    fig.savefig(outpath.replace("upper", tail), bbox_inches="tight", dpi=300)
    logging.info(f"Saved figure {outpath} for {var}_{tail}")


def main(input, output, params):
    EVENTS  = input.events
    FIG1 = output.fig1
    FIG2 = output.fig2
    FIG3 = output.fig3
    FIELDS  = params.fields
    PCRIT   = params.pcrit
    CMAP    = params.cmap
    FIGURES = [FIG1, FIG2, FIG3]

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

    # set values to be floats
    for var in FIELDS.keys():
        for tail_type in ["upper", "lower"]:
            logging.info(f"Statistics for {var}_{tail_type}:")
            logging.info(f"{pk(var, tail_type)}: {gdf[f'{pk(var, tail_type)}'].min()} - {gdf[f'{pk(var, tail_type)}'].max()}")
            logging.info(f"{thresh(var, tail_type)}: {gdf[f'{thresh(var, tail_type)}'].min()} - {gdf[f'{thresh(var, tail_type)}'].max()}")
            logging.info(f"{scale(var, tail_type)}: {gdf[f'{scale(var, tail_type)}'].min()} - {gdf[f'{scale(var, tail_type)}'].max()}")
            logging.info(f"{shape(var, tail_type)}: {gdf[f'{shape(var, tail_type)}'].min()} - {gdf[f'{shape(var, tail_type)}'].max()}")

            logging.info(f"Filling nans for {var}_{tail_type} with -1, 0, 0, 0")
            gdf[pk(var, tail_type)] = gdf[pk(var, tail_type)].fillna(-1)  # switched to -1 to avoid confusion with significant values
            
            logging.info(f"Converting {var}_{tail_type} to float")
            gdf[pk(var, tail_type)] = gdf[pk(var, tail_type)].astype(float)
            gdf[thresh(var, tail_type)] = gdf[thresh(var, tail_type)].astype(float)
            gdf[scale(var, tail_type)] = gdf[scale(var, tail_type)].astype(float)
            gdf[shape(var, tail_type)] = gdf[shape(var, tail_type)].astype(float)

    # turn it into an xarray dataset
    logging.info("Converting GeoDataFrame to xarray Dataset")
    ds = gdf.set_index(['lat', 'lon', 'event']).to_xarray().isel(event=0)

    logging.info(f"Latitude dimensions: {gdf['lat'].nunique()=}")
    logging.info(f"Longitude dimensions: {gdf['lon'].nunique()=}")

    plot_parameters = partial(plot_parameters_for_data, ds=ds,gdf=gdf,
                              fields=FIELDS, cmap=CMAP, pcrit=PCRIT)

    n = 0
    for var in FIELDS.keys():
        for tail in ["upper", "lower"]:
            plot_parameters(var=var, tail=tail, outpath=FIGURES[n])        
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