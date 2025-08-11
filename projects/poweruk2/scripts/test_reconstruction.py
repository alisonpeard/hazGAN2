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
    
    wd = os.path.join("..", "results")
    orig_path = os.path.join(wd, "training", "data.nc")
    ds = xr.open_dataset(orig_path)

    
    # %% make_rgb_images.py & process_generated.py
    if False:
        scaling     = True
        compressing = True
        clipping    = True
        flipping    = False
        save_rgb    = False
        save_npy    = True


        hist_kws = {
            "bins": 50, "density": True, "alpha": 1.,
            "edgecolor": "k", "linewidth": 0.5
        }

        # scaling only:
        # scaling & clipping only: big divergence at {0,1}
        # scaling & compressing & scaling: sorts {0,1} but not divergences around 0.5

        fig, axs = plt.subplots(2, 5, figsize=(20, 8), sharey=True)

        nimgs = [10, ds.time.size][0]
        u_0 = ds.uniform.isel(time=slice(0, nimgs)).values
        axs[0, 0].hist(u_0.ravel(), **hist_kws)
        axs[0, 0].set_title("Original uniform percentiles")

        if flipping:
            u_0 = np.flip(u_0, axis=1)
        
        if clipping:
            u_1 = np.clip(u_0, 1e-6, 1 - 1e-6)

        axs[0, 1].hist(u_1.ravel(), **hist_kws)
        axs[0, 1].set_title("Clipped uniform percentiles")

        y_0 = laplace(u_0)
        axs[0, 2].hist(y_0.ravel(), **hist_kws)
        axs[0, 2].set_title("Laplace percentiles")

        # save the image stats
        array_min = np.min(y_0, axis=(0, 1, 2), keepdims=True)
        array_max = np.max(y_0, axis=(0, 1, 2), keepdims=True)
        n = len(y_0)
        np.savez(os.path.join(tmp_dir, "image_stats.npz"),
                min=array_min, max=array_max, n=n)

        if scaling: # scale to [0, 1]
            y_1 = (y_0 - array_min) / (array_max - array_min)
        else:
            y_1 = y_0
        axs[0, 3].hist(y_1.ravel(), **hist_kws)
        axs[0, 3].set_title("Scaled Laplace percentiles")
        
        if (compressing & scaling): # compress to (0, 1)
            y_2 = (y_1 * (n - 1) + 1) / (n + 1)
        else:
            y_2 = y_1
        axs[0, 4].hist(y_2.ravel(), **hist_kws)
        axs[0, 4].set_title("Compressed Laplace percentiles")

        if save_rgb:
            for i in tqdm(range(nimgs)):
                arr = y_2[i]
                arr = np.uint8(arr * 255)
                img = Image.fromarray(arr, 'RGB')
                output_path = os.path.join(tmp_dir, "rgb", f"footprint{i}.png")
                img.save(output_path)

            # process_generated.py
            _list = glob(os.path.join(tmp_dir, "rgb", "*.png"))
            _images = []
            for png in sorted(_list):
                with Image.open(png) as img:
                    img = np.array(img, dtype=np.float32)
                    _images.append(img)
            y_3 = np.stack(_images, axis=0)
            y_3 /= 255

        elif save_npy:
            for i in tqdm(range(nimgs)):
                arr = y_2[i]
                output_path = os.path.join(tmp_dir, "rgb", f"footprint{i}.npy")
                np.save(output_path, arr)

            _list = glob(os.path.join(tmp_dir, "rgb", "*.npy"))
            _images = []
            for npy in sorted(_list):
                img = np.load(npy)
                _images.append(img)
            y_3 = np.stack(_images, axis=0)

        else:
            y_3 = y_2
        axs[1, 4].hist(y_3.ravel(), **hist_kws)
        axs[1, 4].set_title("Reconstructed Laplace percentiles")

        # apply image statistics to rescale
        image_stats = np.load(os.path.join(tmp_dir, "image_stats.npz"))
        image_minima = image_stats['min']
        image_maxima = image_stats['max']
        image_n      = image_stats['n']

        # rescale images to marginals scale
        if (scaling & compressing):
            y_4 = (y_3 * (image_n + 1) - 1) / (image_n - 1)
        else:
            y_4 = y_3
        axs[1, 3].hist(y_4.ravel(), **hist_kws)
        axs[1, 3].set_title("Decompressed Laplace percentiles")

        if scaling:
            y_4 = y_4 * (image_maxima - image_minima) + image_minima
        axs[1, 2].hist(y_4.ravel(), **hist_kws)
        axs[1, 2].set_title("Re-scaled Laplace percentiles")

        u_2 = inv_laplace(y_4)
        axs[1, 1].hist(u_2.ravel(), **hist_kws)
        axs[1, 1].set_title("Reconstructed uniform percentiles")

        if flipping:
            images_y = np.flip(y_4, axis=1)
            images_u = np.flip(u_2, axis=1)

        axs[1, 0].hist(u_2.ravel(), **hist_kws)
        axs[1, 0].set_title("Reconstructed uniform percentiles (flipped)")



# %% - - - - - - - - dev - - - - - - - -

hist_kws = {
    "bins": 100, "density": True, "alpha": 1.,
    "edgecolor": "k", "linewidth": 0.5
}

if __name__ == "__main__":  # for dev
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
        uniform = np.clip(orig.uniform.values, 1e-6, 1 - 1e-6)  # clip to [0,1]
        orig["laplace"] = (["time", "lat", "lon", "field"], laplace(uniform))

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