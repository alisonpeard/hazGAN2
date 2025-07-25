"""
Examine the daily data extracted by get_data.py.

Loads all the files in the input/ directory.

Plots:
1. Histograms of the overall distributions of wind speed, 30-day rainfall, and wind direction.
2. Sample plots of wind speed, 30-day rainfall, and wind direction for random time steps.
3. Average wind speed, 30-day rainfall, and wind direction over all time.
4. Other statistics (min, max, std) of wind speed and 30-day rainfall.
5. Wind rose plot of wind direction and speed.
"""
# %%
# setup
import os
import glob
import yaml

# for analysis
import numpy as np
import xarray as xr

# for plotting
import cmocean.cm as cmo
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from windrose import WindroseAxes
import matplotlib.pyplot as plt


SUBSET = ["*", "17"][0]

def handle_expver(ds):
    for coord in ["expver", "number"]:
        if coord in ds.dims:
            ds = ds.sel({coord: ds[coord][0]})
        if coord in ds.coords:
            ds = ds.reset_coords(coord, drop=True)
    return ds


def check_missing_days(ds):
    all_days = np.arange(ds.time.min().values, ds.time.max().values + np.timedelta64(1, 'D'), dtype='datetime64[D]')
    missing_days = np.setdiff1d(all_days, ds.time.values)
    if len(missing_days) > 0:
        print(f"Missing days: {len(missing_days)}")
        print(missing_days)
    else:
        print("No missing days found.")


if __name__ == "__main__":

    with open(os.path.join("..", "..", "config.yaml"), "r") as stream:
        config = yaml.safe_load(stream)

    ds_files = glob.glob(os.path.join("..", "..", "results", "processing", "input", f"20{SUBSET}.nc"))
    print(f"Found {len(ds_files)} files in the input directory.")

    ds = xr.open_mfdataset(ds_files)
    ds["from_u"] = - ds.vx * np.sin(ds.dx) # https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
    ds["from_v"] = - ds.vx * np.cos(ds.dx)
    ds["to_u"] = - ds.from_u
    ds["to_v"] = - ds.from_v
    ds = ds.sortby(["latitude", "longitude"])
    check_missing_days(ds)

    # %% histograms of overall distributions
    fig, axs = plt.subplots(1, 3, figsize=(10, 5))

    hist_kws = {"color": "lightgrey", "linewidth": 0.5, "edgecolor": "k",
                "bins": 50, "density": True}

    ds.vx.plot.hist(ax=axs[0], **hist_kws)
    axs[0].set_xlabel("Wind speed (vx) [m/s]")

    ds.r30.plot.hist(ax=axs[1], **hist_kws)
    axs[1].set_xlabel("30-day rainfall (r3) [m]")

    ds.dx.plot.hist(ax=axs[2], **hist_kws)
    axs[2].set_xlabel("Wind direction (dx) [degrees]")

    plt.tight_layout()

    # %% plot some samples
    nrows = 5

    fig, axes = plt.subplots(nrows, 2, subplot_kw={"projection": ccrs.PlateCarree()},
                                figsize=(16, 6 * nrows))
    for row in range(nrows):
        axs = axes[row, :]
        t = np.random.randint(0, ds.sizes["time"])

        ds.isel(time=t).vx.plot(ax=axs[0], cmap=cmo.speed, cbar_kwargs={"label": "max gust speed (m/s)", "shrink": 0.6})
        ds.isel(time=t).r30.plot(ax=axs[1], cmap=cmo.rain, cbar_kwargs={"label": "30-day rainfall (m)", "shrink": 0.6})

        resample = ds.isel(time=t, longitude=slice(None, None, 7), latitude=slice(None, None, 7))
        resample.plot.streamplot(x='longitude', y='latitude', u='to_u', v='to_v', 
                                    transform=ccrs.PlateCarree(), color="k", ax=axs[0], density=1.5,
                                    linewidth=0.5, arrowstyle='->', arrowsize=1.25
                                    )
        for ax in axs:
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=1)
            ax.set_extent([
                ds.longitude.min().values,
                ds.longitude.max().values,
                ds.latitude.min().values,
                ds.latitude.max().values
            ])

        time_str = ds.time[t].values.astype("datetime64[s]").astype(str)
        axs[0].set_title(f"ERA5 wind speed and direction (t={time_str})")
        axs[1].set_title(f"ERA5 30-day rainfall (t={time_str})")

    fig.suptitle("ERA5 data samples", fontsize=16, y=1.0)
    plt.tight_layout()
    plt.show()

    # %% visualise average over all time
    ds_mean = ds.mean(dim="time")
    fig, axs = plt.subplots(1, 2, subplot_kw={"projection": ccrs.PlateCarree()},
                            figsize=(14, 6), sharex=True, sharey=True)
    ds_mean.vx.plot.contourf(ax=axs[0], cmap=cmo.speed, levels=50, cbar_kwargs={"label": "max gust speed (m/s)", "shrink": 0.6})
    ds_mean.r30.plot.contourf(ax=axs[1], cmap=cmo.rain, levels=50, cbar_kwargs={"label": "30-day rainfall (m)", "shrink": 0.6})
    resample = ds_mean.isel(longitude=slice(None, None, 7), latitude=slice(None, None, 7))
    streamlines = resample.plot.streamplot(x='longitude', y='latitude', u='to_u', v='to_v', 
                                transform=ccrs.PlateCarree(), color="k", ax=axs[0], density=1.5,
                                linewidth=0.5, arrowstyle='->', arrowsize=1.25
                                )
    for ax in axs:
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.75)
        ax.set_extent([
            ds_mean.longitude.min().values,
            ds_mean.longitude.max().values,
            ds_mean.latitude.min().values,
            ds_mean.latitude.max().values
        ])

    axs[0].set_title("ERA5 wind speed and direction")
    axs[1].set_title("ERA5 30-day rainfall")
    fig.suptitle("ERA5 data averaged over all time", fontsize=16)
    plt.tight_layout()
    plt.show()

    # %% visualise other statistics
    for stat in ["min", "max", "std"]:
        ds_stat = getattr(ds, stat)(dim="time")
        fig, axs = plt.subplots(1, 2, subplot_kw={"projection": ccrs.PlateCarree()},
                                figsize=(14, 6), sharex=True, sharey=True)
        ds_stat.vx.plot.contourf(ax=axs[0], cmap=cmo.speed, levels=50, cbar_kwargs={"label": f"max gust speed ({stat}) (m/s)", "shrink": 0.6})
        ds_stat.r30.plot.contourf(ax=axs[1], cmap=cmo.rain, levels=50, cbar_kwargs={"label": f"30-day rainfall ({stat}) (m)", "shrink": 0.6})
        resample = ds_stat.isel(longitude=slice(None, None, 7), latitude=slice(None, None, 7))

        for ax in axs:
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.75)
            ax.set_extent([
                ds_stat.longitude.min().values,
                ds_stat.longitude.max().values,
                ds_stat.latitude.min().values,
                ds_stat.latitude.max().values
            ])

        axs[0].set_title(f"ERA5 wind speed and direction ({stat})")
        axs[1].set_title(f"ERA5 30-day rainfall ({stat})")
        fig.suptitle(f"ERA5 data {stat} over all time", fontsize=16)
        plt.tight_layout()
        plt.show()

    # %% wind rose overall etc.
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='windrose'))

    ds_agg = ds.mean(dim="time")

    dirs = ds_agg.dx.values.flatten()
    speeds = ds_agg.vx.values.flatten()

    ax.bar(dirs, speeds, normed=True, opening=0.8, edgecolor='black',
           nsector=16, bins=10,
           cmap=cmo.speed, linewidth=0.5)
    ax.set_legend()
    
# %%  

