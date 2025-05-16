"""Load all generated PNGs, transform to original scales, and save to NetCDF."""
# %%
import os
from glob import glob
import logging
from snakemake.script import snakemake

os.environ["USE_PYGEOS"] = "0"
from PIL import Image
import numpy as np
import xarray as xr

# os.system("python -m pip install ../../packages/hazGAN/")

from hazGAN.statistics import gumbel, inv_gumbel, invPIT

if __name__ == "__main__":
    # configure logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    
    # load parameters
    IMAGE_DIR   = snakemake.input.image_dir
    IMAGE_STATS = snakemake.input.image_stats
    TRAIN       = snakemake.input.training_data
    DO_SUBSET   = snakemake.params.do_subset # NOTE: this is messy... clean later
    THRESH      = snakemake.params.event_subset
    FIELDS      = snakemake.params.fields
    OUTPUT      = snakemake.output.netcdf

    # load all images
    IMAGES = glob(os.path.join(IMAGE_DIR, "*.png"))
    logging.info(f"Found {len(IMAGES)} images in {IMAGE_DIR}.")
    images = []
    for png in sorted(IMAGES):
        img = Image.open(png)
        images.append(img)
    images = np.array(images, dtype=np.float32)
    images /= 255
    images = np.flip(images, axis=1) # flip y-axis (latitude)
    logging.info(f"Created generated images ndarray of shape {images.shape}.")

    # apply image statistics to rescale
    image_stats = np.load(IMAGE_STATS)
    image_minima = image_stats['min']
    image_maxima = image_stats['max']
    image_n      = image_stats['n']
    image_range  = image_maxima - image_minima
    logging.info(f"Image statistics: min {image_minima}, max {image_maxima}, n {image_n}.")

    # rescale images to Gumbel scale
    images_g = (images * (image_n + 1) - 1) / (image_n - 1) * image_range + image_minima
    images_u = inv_gumbel(images_g) # np.exp(-np.exp(-images_g))

    # load the training data for the inverse ECDF and parameters
    # NOTE: here is a good place to check distribution of data... skipping for now
    train = xr.open_dataset(TRAIN)
    
    # subset by threshold if needed
    train['intensity'] = getattr(train.sel(field=THRESH["field"]).anomaly, THRESH["func"])(dim=['lon', 'lat'])
    if DO_SUBSET:
        mask = (train['intensity'] > THRESH["value"]).values
        idx  = np.where(mask)[0]
        train   = train.isel(time=idx)

    train   = train.sortby("intensity", ascending=False)
    train_x = train["anomaly"].values
    train_u = train["uniform"].values
    train_g = gumbel(train_u)
    params  = train["params"].values
    logging.info(f"Loaded anomaly training data of shape {train_x.shape}.")
    logging.info(f"Loaded uniform traingin data of shape {train_u.shape}.")
    logging.info(f"Loaded parameters of shape {params.shape}.")

    # check range of uniform samples
    if not ((train_u.max() < 1.) and (train_u.min() > 0.)):
        logging.error("Training uniform samples not in [0,1] range")
        logging.error(f"Max: {train_u.max()}, Min: {train_u.min()}")
        train_u = np.clip(train_u, 1e-6, 1 - 1e-6)
    
    if not ((images_u.max() < 1.) and (images_u.min() > 0.)):
        logging.error("Generated uniform samples not in [0,1] range")
        logging.error(f"Max: {images_u.max()}, Min: {images_u.min()}")
        images_u = np.clip(images_u, 1e-6, 1 - 1e-6)

    # transform images to original scale using invPIT
    # TODO: need to check y-axis orientation
    distns = [field['distn'] for field in FIELDS.values()]
    images_x = invPIT(images_u, train_x, params, distns=distns)

    # save to NetCDF
    data = xr.Dataset(
        {
            "anomaly": (("time", "lat", "lon", "field"), images_x),
            "uniform": (("time", "lat", "lon", "field"), images_u),
            "params": (("lat", "lon", "param", "field"), params),
        },
        coords={
            "param": ["loc", "scale", "shape"],
            "field": list(FIELDS.keys()),
            "time": (("time"), np.arange(images_x.shape[0])),
            "lat": (("lat"), train["lat"].values),
            "lon": (("lon"), train["lon"].values),
        },
    )

    data.to_netcdf(OUTPUT, format="NETCDF4")
    logging.info(f"Saved generated data to {OUTPUT}.")
    logging.info("Done.")
# %%
