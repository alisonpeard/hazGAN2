"""Check reconstructed train matches original."""
# %%
import os
import yaml
import shutil

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from PIL import Image
from glob import glob
from warnings import warn
from tqdm import tqdm


def load_images(path):
    image_list = glob(os.path.join(path, "*.png"))
    print(f"Found {len(image_list)} images in {path}.")
    images = []
    for png in sorted(image_list):
        with Image.open(png) as img:
            img = np.array(img, dtype=np.float32)
            images.append(img)
    images = np.stack(images, axis=0)
    images /= 255
    print(f"Created images ndarray of shape {images.shape}.")
    return images


def laplace(uniform, mu=0, b=1):
    """uniform -> Laplace(mu, b) (quantile function)"""
    maxval = np.max(uniform)
    if maxval == 1:
        warn("Values == 1 found, scaling by 1e-6")
        uniform *= 1 - 1e-6
    if maxval > 1:
        raise ValueError(f"Some uniform > 1 ({maxval})")
    
    return np.where(
        uniform <= 0.5, 
        mu + b * np.log(2 * uniform),
        mu - b * np.log(2 - 2 * uniform)
        )


def inv_laplace(x, mu=0, b=1):
    """Laplace(mu, b) -> uniform (CDF function)."""
    return np.where(
        x <= mu,
        0.5 * np.exp((x - mu) / b),
        1 - 0.5 * np.exp(-(x - mu) / b)
    )


hist_kws = {
    "bins": 100, "density": True, "alpha": 0.5,
    "edgecolor": "k", "linewidth": 0.5
    }


if __name__ == "__main__":
    # mfe for making images and reconstructing
    tmp_dir = os.path.join("..", "tmp")

    # delete old tmp dir if exists
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    os.makedirs(os.path.join(tmp_dir, "rgb"))
    
    wd = os.path.join("..", "results")
    orig_path = os.path.join(wd, "training", "data.nc")
    ds = xr.open_dataset(orig_path)

    
    # %% make_rgb_images.py
    scaling     = True
    compressing = False
    clipping    = True
    flipping    = False
    save_rgb    = False

    # scaling only:
    # scaling & clipping only: big divergence at {0,1}
    # scaling & compressing & scaling: sorts {0,1} but not divergences around 0.5

    nimgs = [10, ds.time.size][0]
    array_u = ds.uniform.isel(time=slice(0, nimgs)).values

    if flipping:
        array_u = np.flip(array_u, axis=1)
    
    if clipping:
        array_u = np.clip(array_u, 1e-6, 1 - 1e-6)

    array_y = laplace(array_u)

    # save the image stats
    array_min = np.min(array_y, axis=(0, 1, 2), keepdims=True)
    array_max = np.max(array_y, axis=(0, 1, 2), keepdims=True)
    n = len(array_y)
    np.savez(os.path.join(tmp_dir, "image_stats.npz"),
             min=array_min, max=array_max, n=n)

    if scaling: # scale to [0, 1]
        array_y = (array_y - array_min) / (array_max - array_min)
    else:
        array_y = array_y
    
    if (compressing & scaling): # compress to (0, 1)
        array_y = (array_y * (n - 1) + 1) / (n + 1)
    else:
        array_y = array_y

    if save_rgb:
        for i in tqdm(range(nimgs)):
            arr = array_y[i]
            arr = np.uint8(arr * 255)
            img = Image.fromarray(arr, 'RGB')
            output_path = os.path.join(tmp_dir, "rgb", f"footprint{i}.png")
            img.save(output_path)

        # process_generated.py
        train_list = glob(os.path.join(tmp_dir, "rgb", "*.png"))
        images_train = []
        for png in sorted(train_list):
            with Image.open(png) as img:
                img = np.array(img, dtype=np.float32)
                images_train.append(img)
        images_train = np.stack(images_train, axis=0)
        images_train /= 255
    else:
        images_train = array_y

    # apply image statistics to rescale
    image_stats = np.load(os.path.join(tmp_dir, "image_stats.npz"))
    image_minima = image_stats['min']
    image_maxima = image_stats['max']
    image_n      = image_stats['n']

    # rescale images to marginals scale
    if (scaling & compressing):
        images_y = (images_train * (image_n + 1) - 1) / (image_n - 1)
    else:
        images_y = images_train

    if scaling:
        images_y = images_y * (image_maxima - image_minima) + image_minima

    inv_transform = inv_laplace
    images_u = inv_transform(images_y)

    if flipping:
        images_y = np.flip(images_y, axis=1)
        images_u = np.flip(images_u, axis=1)

    # %compare the original and reconstructed uniform distributions
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(array_u.flatten(), label='Training', **hist_kws)
    ax.hist(images_u.flatten(), label='Reconstructed', **hist_kws)
    ax.set_title("Histogram of percentiles")
    ax.set_xlabel("Percentile")
    ax.set_ylabel("Density")
    ax.legend()
    plt.show()

    # %% compare the original and reconstructed uniform distributions
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(array_y.flatten(), label='Training', **hist_kws)
    ax.hist(images_y.flatten(), label='Reconstructed', **hist_kws)
    ax.set_title("Histogram of percentiles")
    ax.set_xlabel("Percentile")
    ax.set_ylabel("Density")
    ax.legend()
    plt.show()

    # %%
    # Simulate the 8-bit conversion
    quantized = np.uint8(array_y * 255) / 255.0
    quantization_error = np.mean(np.abs(array_y - quantized))
    print(f"Quantization error: {quantization_error}")


# %% - - - - - - - - dev - - - - - - - -
if __name__ == "__main__" and False:  # for dev
    # load everything
    wd = os.path.join("..", "results")
    orig_path = os.path.join(wd, "training", "data.nc")
    recon_path = os.path.join(wd, "generated", "netcdf", "train.nc")
    image_path = os.path.join(wd, "training", "rgb")
    stats_path = os.path.join(wd, "training", "image_stats.npz")

    cfg_path = os.path.join("..", "config.yaml")
    with open(cfg_path, "r") as stream:
        cfg = yaml.safe_load(stream)

    orig = xr.open_dataset(orig_path)
    recon = xr.open_dataset(recon_path)

    images = load_images(image_path)
    image_stats = np.load(stats_path)
    image_min = image_stats['min']
    image_max = image_stats['max']
    image_n = image_stats['n']

    #%% laplace transforms
    uniform = np.clip(orig.uniform.values, 1e-6, 1 - 1e-6)  # clip to [0,1]
    orig["laplace"] = (["time", "lat", "lon", "field"],
                  laplace(uniform))

    # %%
    # check: can I recover orig::uniform from images
    images_y = (images * (image_n + 1) - 1) / (image_n - 1)
    images_y = images_y * (image_max - image_min) + image_min
    images_u = inv_laplace(images_y)

    # %%
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(orig["laplace"].values.flatten(), label='Training', **hist_kws)
    ax.hist(images_y.flatten(), label='Reconstructed', **hist_kws)
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(orig["uniform"].values.flatten(), label='Training', **hist_kws)
    ax.hist(images_u.flatten(), label='Reconstructed', **hist_kws)
    plt.show()


    # %%
    domain = "uniform"

    for field, properties in cfg["fields"].items():
        title = properties.get("title", field)
        x = orig[domain].sel(field=field).values.flatten()
        x_recon = recon[domain].sel(field=field).values.flatten()

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(x, label='Training', **hist_kws)
        ax.hist(x_recon, label='Reconstructed', **hist_kws)
        ax.set_title(f"Histogram of {title.lower()} percentiles")
        ax.set_xlabel(f"{title} anomaly")
        ax.set_ylabel("Density")
        ax.legend()
        plt.show()

    # %%
    domain = "laplace"

    hist_kws = {
        "bins": 100, "density": True, "alpha": 0.5,
        "edgecolor": "k", "linewidth": 0.5
        }

    for field, properties in cfg["fields"].items():
        title = properties.get("title", field)
        x = orig[domain].sel(field=field).values.flatten()
        x_recon = recon[domain].sel(field=field).values.flatten()

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(x, label='Training', **hist_kws)
        ax.hist(x_recon, label='Reconstructed', **hist_kws)
        ax.set_title(f"Histogram of {title.lower()} percentiles")
        ax.set_xlabel(f"{title} anomaly")
        ax.set_ylabel("Density")
        ax.legend()
        plt.show()

    # %%
    domain = "anomaly"

    hist_kws = {
        "bins": 100, "density": True, "alpha": 0.5,
        "edgecolor": "k", "linewidth": 0.5
        }

    for field, properties in cfg["fields"].items():
        title = properties.get("title", field)
        x = orig[domain].sel(field=field).values.flatten()
        x_recon = recon[domain].sel(field=field).values.flatten()

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(x, label='Training', **hist_kws)
        ax.hist(x_recon, label='Reconstructed', **hist_kws)
        ax.set_title(f"Histogram of {title.lower()} anomalies")
        ax.set_xlabel(f"{title} anomaly")
        ax.set_ylabel("Density")
        ax.legend()
        plt.show()

    # %%
    # forward transform
    array_u = orig.uniform.values
    array_u = np.flip(array_u, axis=1) # flip latitude
    array_u = np.clip(array_u, 1e-6, 1 - 1e-6) # clip to [0,1]
    array_y = laplace(array_u)

    array_min = np.min(array_y, axis=(0, 1, 2), keepdims=True)
    array_max = np.max(array_y, axis=(0, 1, 2), keepdims=True)
    n = len(array_y)

    array = (array_y - array_min) / (array_max - array_min)
    array = (array * (n - 1) + 1) / (n + 1)

    # reverse transform
    images_y = (array * (n + 1) - 1) / (n - 1)
    images_y = images_y * (array_max - array_min) + array_min
    images_y = np.flip(images_y, axis=1)

    images_u = inv_laplace(images_y)
    images_u = np.flip(images_u, axis=1)

    # two eval plots
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(array_y.flatten(), bins=100, density=True, alpha=0.5,
            edgecolor="k", linewidth=0.5, label='Training')
    ax.hist(images_y.flatten(), bins=100, density=True, alpha=0.5,
            edgecolor="k", linewidth=0.5, label='Reconstructed')
    ax.set_title("Histogram of Laplace percentiles")
    ax.set_xlabel("Laplace anomaly")
    ax.set_ylabel("Density")
    ax.legend()
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(array_u.flatten(), bins=100, density=True, alpha=0.5,
            edgecolor="k", linewidth=0.5, label='Training')
    ax.hist(images_u.flatten(), bins=100, density=True, alpha=0.5,
            edgecolor="k", linewidth=0.5, label='Reconstructed')
    ax.set_title("Histogram of percentiles")
    ax.set_xlabel("Percentile")
    ax.set_ylabel("Density")
    ax.legend()
    plt.show()

# %%