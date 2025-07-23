"""Test the get_data.py results for the year 2020."""
# %%
import os
import yaml
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

if __name__ == "__main__":

    with open(os.path.join("..", "..", "config.yaml"), "r") as stream:
        config = yaml.safe_load(stream)

    ds_path = os.path.join("..", "..", "results", "processing", "input", "2017.nc")
    ds = xr.open_dataset(ds_path)

    # for plotting
    ds["u"] = - ds.vx * np.sin(ds.dx)
    ds["v"] = - ds.vx * np.cos(ds.dx)

    ds = ds.sortby(["latitude", "longitude"])

    # visualise the data
    for _ in range(10):
        t = np.random.randint(0, ds.sizes["time"])

        fig, axs = plt.subplots(1, 2, subplot_kw={"projection": ccrs.PlateCarree()},
                                figsize=(20, 8))
        
        ds.isel(time=t).vx.plot(ax=axs[0], cmap="Spectral_r", cbar_kwargs={"label": "max gust speed (m/s)", "shrink":0.8})
        ds.isel(time=t).r30.plot(ax=axs[1], cmap="PuBu", cbar_kwargs={"label": "30-day rainfall (m)", "shrink":0.8})

        resample = ds.isel(time=t,longitude=slice(None, None, 7), latitude=slice(None, None, 7))
        streamlines = resample.plot.streamplot(x='longitude', y='latitude', u='u', v='v', 
                                    transform=ccrs.PlateCarree(), color="k", ax=axs[0], density=1.5,
                                    linewidth=0.5, arrowstyle='->', arrowsize=0.5
                                    )
        for ax in axs:
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=1)

        axs[0].set_title("ERA5 wind speed and direction")
        axs[1].set_title("ERA5 30-day rainfall")
        fig.suptitle(f"ERA5 data for {ds.time[t].values}", fontsize=16)
        plt.show()

# %%  

