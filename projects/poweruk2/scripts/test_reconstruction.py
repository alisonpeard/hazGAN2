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


def load_numpy_images(path, n=10):
    _list = glob(os.path.join(path, "*.npy"))
    print(f"Found {len(_list)} numpy images in {path}.")
    images = []
    for f in tqdm(sorted(_list[:n])):
        img = np.load(f)
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


if __name__ == "__main__":
    # mfe for making images and reconstructing
    tmp_dir = os.path.join("..", "tmp")

    # delete old tmp dir if exists
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    os.makedirs(os.path.join(tmp_dir, "rgb"))


    # %% - - - - - - - - dev - - - - - - - -

    hist_kws = {
        "bins": 100, "density": True, "alpha": 1.,
        "edgecolor": "k", "linewidth": 0.5
    }

    if True:
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

        images = load_numpy_images(image_path)
        image_stats = np.load(stats_path)
        image_min = image_stats['min']
        image_max = image_stats['max']
        image_n = image_stats['n']

        # %% laplace transforms
        uniform = np.clip(orig["uniform"].values, 1e-6, 1 - 1e-6)  # clip to [0,1]
        orig["standardised"] = (["time", "lat", "lon", "field"], laplace(uniform))

        # %%
        # check: can I recover orig::uniform from images
        images_y = (images * (image_n + 1) - 1) / (image_n - 1)
        images_y = images_y * (image_max - image_min) + image_min
        images_u = inv_laplace(images_y)

        # # %%
        # fig, ax = plt.subplots(figsize=(10, 5))
        # ax.hist(orig["laplace"].values.flatten(), label='Training', **hist_kws)
        # ax.hist(images_y.flatten(), label='Reconstructed', **hist_kws)
        # plt.show()

        # fig, ax = plt.subplots(figsize=(10, 5))
        # ax.hist(orig["uniform"].values.flatten(), label='Training', **hist_kws)
        # ax.hist(images_u.flatten(), label='Reconstructed', **hist_kws)
        # plt.show()

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
        domain = "standardised"

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

        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

# %%