"""
Rain bombs are a known issue in ERA5 data[2]. Though not addressed in the documentation,
the Bay of Bengal wind data has clusters of very high winds that do not look realistic[1].
This script detects these clusters using blob detection and visualizes them.

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
- 1997-05-16
- 1992-04-15
- 1952-05-09

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

THRESH = 0.015
SIGMA = 2
INPUT_DIR = "/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/"
PROBLEM_DATES = [
    "1952-05-09",
    "1992-04-15",
    "1997-05-16"
]
SAMPLES = ['/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1997.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1992.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1952.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1995.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1963.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1953.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1950.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1949.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1946.nc',
                '/Volumes/mistral/alison/hazGAN2/projects/bayofbengal_era5/results/processing/input/1945.nc'
]

def has_blobs(x, max_sigma=SIGMA, threshold=THRESH):
    """Apply blob detection to the dataset."""
    blobs = blob_doh(x, max_sigma=max_sigma, threshold=threshold)
    return len(blobs)


def get_blobs(x, max_sigma=SIGMA, threshold=THRESH):
    """Apply blob detection to the dataset."""
    blobs = blob_doh(x, max_sigma=max_sigma, threshold=threshold)
    return blobs


if __name__ == "__main__":
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
        ax.imshow(arr)
        for blob in blobs:
            y, x, r = blob
            c = plt.Circle((x, y), r, color="lime", linewidth=2, fill=False)
            ax.add_patch(c)
        ax.set_axis_off()
        plt.show()


    # %%
    files = sorted(glob(os.path.join(INPUT_DIR, "*.nc")), reverse=True)
    print(f"Found {len(files)} files in {INPUT_DIR}.")

    filtered = []
    for file in (pbar := tqdm(files, desc="Processing files")):
        year = os.path.basename(file).split('.')[0]
        pbar.set_postfix({"year": year})
        try:
            ds = xr.open_dataset(file)
            ds["ws_norm"] = (ds["ws"] - ds["ws"].min()) / (ds["ws"].max() - ds["ws"].min())
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
