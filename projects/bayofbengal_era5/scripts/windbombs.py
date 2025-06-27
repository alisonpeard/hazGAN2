"""
Rain bombs are a known issue in ERA5 data[2]. Though not addressed in the documentation,
the Bay of Bengal wind data has clusters of very high winds that do not look realistic[1].
This script detects these clusters using blob detection and visualizes them.

This is messy, hard to get the right parameters so using these to manually choose exclusion dates.

Last run (2025-06-24) detected blobs in the following dates:
- 2002-05-10
- 1997-05-16
- 1992-04-15
- 1984-06-21
- 1952-05-09
- 1946-05-12
- 1946-05-13
- 1946-06-06
- 1946-06-07
- 1946-07-20
- 1946-09-08
- 1945-06-07

From manual inspection, the following dates were chosen to exclude:
  - "1952-05-09"
  - "1992-04-15"
  - "1997-05-16"
  - "2008-05-04"
  - "1992-05-25"
  - "1997-05-24"
  - "1995-05-02"

NOTE: 1940 was omitted as it was corrupted on the SoGE cluster.

References: 
..[1] https://forum.ecmwf.int/t/unusual-10m-wind-distribution-over-bay-of-bengal/7399 
..[2] https://forum.ecmwf.int/t/era5-versus-era-interim-total-precipitation/978
 
"""
# %%
import os
import pandas as pd
from glob import glob
from tqdm import tqdm
import xarray as xr
import matplotlib.pyplot as plt
from skimage.feature import blob_doh

THRESH = 0.01 #0.015
SIGMA = 3
INPUT_DIR = "/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/"
PROBLEM_DATES = [
    # original
    "1952-05-09",
    "1992-04-15",
    "1997-05-16",
    # found after training
    "2008-05-04", # only detected when THRESH=0.01 and SIGMA=3
    "1992-05-25", # needs THRESH<=0.00125
    "1997-05-24", # needs THRESH<=0.0075
    "1995-05-02", # THRESH=0.0025
]
SAMPLES = [
    # 2020,
    # 1997,
    # 1952,
    1992,
    1995,
    1997,
    2008
]

SAMPLES = [os.path.join(INPUT_DIR, f"{year}.nc") for year in SAMPLES]
# SAMPLES = ["/Users/alison/Downloads/1992.nc"]

def has_blobs(x, max_sigma=SIGMA, threshold=THRESH):
    """Apply blob detection to the dataset."""
    blobs = blob_doh(x, max_sigma=max_sigma, threshold=threshold)
    return len(blobs)


def get_blobs(x, max_sigma=SIGMA, threshold=THRESH):
    """Apply blob detection to the dataset."""
    blobs = blob_doh(x, max_sigma=max_sigma, threshold=threshold)
    return blobs

# %%
if __name__ == "__main__":
    wdir = os.path.dirname(__file__)
    # run this first to manually check the parameters are finding the right blobs
    for file in tqdm(SAMPLES, desc="Processing sample files"):
        ds = xr.open_dataset(file)
        print(file)
        ds["ws_norm"] = (ds["ws"] - ds["ws"].min()) / (ds["ws"].max() - ds["ws"].min())

        for t in range(len(ds.time)):
            arr = ds.isel(time=t).ws_norm.values
            blobs = get_blobs(arr)
            if len(blobs) > 0:
                break


        print(f"Found {len(blobs)} for {file} at time {ds.time[t].values}.")
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.imshow(arr, cmap="binary_r")
        for blob in blobs:
            y, x, r = blob
            c = plt.Circle((x, y), r, color="r", linewidth=1, fill=False)
            ax.add_patch(c)
        ax.set_axis_off()
        plt.show()

        fig.savefig(os.path.join("..", "results", "figures", f"windbomb_{t}.png"), dpi=300, bbox_inches="tight")

    # %% manually inspect known problem dates to calibrate parameters
    for date in PROBLEM_DATES:
        # idx = 4
        # date = PROBLEM_DATES[-idx]
        year = date.split("-")[0]
        file = os.path.join(INPUT_DIR, f"{year}.nc")
        ds = xr.open_dataset(file)
        print(f"Processing {file} for date {date}.")
        
        # normalise each day independently
        ds["ws_min"] = ds["ws"].min(dim=["longitude", "latitude"])
        ds["ws_max"] = ds["ws"].max(dim=["longitude", "latitude"])
        ds["ws_norm"] = (ds["ws"] - ds["ws_min"]) / (ds["ws_max"] - ds["ws_min"])
        ds.sel(time=date).ws.plot(figsize=(3, 3), cmap="binary_r")
        plt.show()

        arr = ds.sel(time=date).ws_norm.values
        blobs = get_blobs(arr, max_sigma=3, threshold=0.01)
        if len(blobs) > 0:
            print(f"Found {len(blobs)} blobs for {date}.")
            fig, ax = plt.subplots(figsize=(3, 3))
            ax.imshow(arr, cmap="binary_r")
            for blob in blobs:
                y, x, r = blob
                c = plt.Circle((x, y), r, color="r", linewidth=1, fill=False)
                ax.add_patch(c)
            ax.set_axis_off()
            plt.show()
            fig.savefig(os.path.join("..", "results", "figures", f"windbomb_{date}.png"), dpi=300, bbox_inches="tight")

        

        # %%
        if False: # debugging
            data = xr.open_dataset("/Users/alison/Desktop/data.nc")
            import numpy as np
            metadata = pd.read_parquet("/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/event_metadata.parquet")
            metadata["time"] = pd.to_datetime(metadata["time"])
            metadata["year"] = metadata["time"].dt.year
            metadata = metadata[metadata["time"] == pd.to_datetime(date)].sort_values("time")
            start = metadata["time"].min()
            end = metadata["time"].max()

            data.sel(time=date, field="ws").anomaly.plot()
            date_list = pd.date_range(start, end, freq="D")
            fig, axs = plt.subplots(1, len(date_list), figsize=(3 * len(date_list), 3))
            axs = axs.flatten() if len(date_list) > 1 else [axs]
            for i, date_i in enumerate(date_list):
                ds.sel(time=date_i).ws.plot(ax=axs[i])

            ds_footprint = ds.sel(time=slice(start, end))
            ds_footprint["ws"] = ds_footprint.ws.max(dim="time")
            ds_footprint.ws.plot()

            #  now argmax of maximum wind speed
            ds_argmax = ds.sel(time=slice(start, end))
            ds_argmax["ws"] = ds_argmax.ws.argmax(dim="time").astype(int)

            vmin = ds_argmax.ws.min().values
            vmax = ds_argmax.ws.max().values
            levels = np.arange(vmin, vmax + 2)
            ds_argmax.ws.plot(levels=levels)

            # ds.sel(time=date_list[-1]).ws.plot()
            ds.sel(time=date).ws_norm.plot()
            arr = ds.sel(time=date).ws_norm.values

            sigma = 3 #SIGMA
            threshold = 0.0025 #THRESH

            blobs = get_blobs(arr, max_sigma=sigma, threshold=threshold)
            print(f"Found {len(blobs)} blobs for {date} with {sigma=} and {threshold=}.")
            if len(blobs) > 0:
                print(f"Found {len(blobs)} blobs for {date}.")
                fig, ax = plt.subplots(figsize=(3, 3))
                ax.imshow(arr, cmap="binary_r")
                for blob in blobs:
                    y, x, r = blob
                    c = plt.Circle((x, y), r, color="r", linewidth=1, fill=False)
                    ax.add_patch(c)
                ax.set_axis_off()
                plt.show()


    # %% then run this to process all files
    files = sorted(glob(os.path.join(INPUT_DIR, "*.nc")), reverse=True)
    print(f"Found {len(files)} files in {INPUT_DIR}.")

    filtered = []
    for file in (pbar := tqdm(files, desc="Processing files")):
        year = os.path.basename(file).split('.')[0]
        pbar.set_postfix({"year": year})
        try:
            ds = xr.open_dataset(file)
            # normalise each day independently
            ds["ws_min"] = ds["ws"].min(dim=["longitude", "latitude"])
            ds["ws_max"] = ds["ws"].max(dim=["longitude", "latitude"])
            ds["ws_norm"] = (ds["ws"] - ds["ws_min"]) / (ds["ws_max"] - ds["ws_min"])
            blobs = xr.apply_ufunc(
                has_blobs, 
                ds.ws_norm, 
                input_core_dims=[['latitude', 'longitude']],  # consume spatial dims
                output_core_dims=[[]],             # return scalar
                vectorize=True
            )
            mask = blobs > 0
            if sum(mask) > 0:
                filtered_ds = ds.isel(time=mask)
                filtered.append(filtered_ds)
                print(f"Found {blobs.sum().data} blobs.")
            else:
                print(f"Found no blobs.")
        except Exception as e:
            print(f"Error processing {file}: {e}")
            continue

    # %%
    filtered_all = xr.concat(filtered, dim="time")

    for t in range(len(filtered_all.time)):
        fig, ax = plt.subplots()
        plotme = filtered_all.isel(time=t)
        plotme.ws.plot(ax=ax)
        plt.show()

    dates = filtered_all.time.values
    dates = [pd.to_datetime(date).strftime("%Y-%m-%d") for date in dates]
    for date in dates:
        print(date)


# %%
