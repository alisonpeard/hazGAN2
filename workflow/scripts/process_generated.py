"""Load all generated PNGs, inverse-PIT to original scales, and save to NetCDF."""
import os
import logging
from glob import glob

os.environ["USE_PYGEOS"] = "0"

from PIL import Image
import numpy as np
import xarray as xr

from src import funcs
import src.python.statistics as stats


def subset_func(ds:xr.Dataset, subset:dict):
    """Subset the dataset using function and threshold."""
    func = getattr(funcs, subset["func"])
    args = subset["args"]
    thresh = subset["value"]

    logging.info(f"Subsetting {func}{*args,} with threshold {thresh}.")

    for arg in args:
        ds[arg] = ds.sel(field=arg).anomaly

    intensity = func(ds, *args, dim=["lon", "lat"])

    for arg in args:
        ds = ds.drop_vars(arg)

    mask = (intensity > thresh).values
    idx = np.where(mask)[0]
    return ds.isel(time=idx)


def main(input, output, params):
    # load all images
    flist = glob(os.path.join(input.image_dir, "*.npy"))
    logging.info(f"Found {len(flist)} images in {input.image_dir}.")
    images = []
    for f in sorted(flist):
        img = np.load(f)
        images.append(img / 255)
    images = np.stack(images, axis=0)
    logging.info(f"Created generated images ndarray of shape {images.shape}.")

    # load training images
    flist = glob(os.path.join(input.training_dir, "*.npy"))
    logging.info(f"Found {len(flist)} training images in {input.training_dir}.")
    images_train = []
    for f in sorted(flist):
        img = np.load(f)
        images_train.append(img / 255)
    images_train = np.stack(images_train, axis=0)
    logging.info(f"Created training images ndarray of shape {images_train.shape}.")

    # apply image statistics to rescale
    image_stats = np.load(input.image_stats)
    image_minima = image_stats['min']
    image_maxima = image_stats['max']
    # image_n      = image_stats['n']
    image_range  = image_maxima - image_minima
    logging.info(f"Image statistics: min {image_minima}, max {image_maxima}") #, n {image_n}.")

    inv_transform = getattr(stats, f"inv_{params.domain}")

    # rescale images to marginals scale
    # images_y = (images * (image_n + 1) - 1) / (image_n - 1)

    images_y = images * image_range + image_minima
    images_u = inv_transform(images_y)
    images_y = np.flip(images_y, axis=1)
    images_u = np.flip(images_u, axis=1)

    # rescale train in same way
    # compare_y = (images_train * (image_n + 1) - 1) / (image_n - 1)
    compare_y = images_train * image_range + image_minima
    compare_u = inv_transform(compare_y)
    compare_y = np.flip(compare_y, axis=1)
    compare_u = np.flip(compare_u, axis=1)

    # load training netCDF
    train = xr.open_dataset(input.training_data)
    
    if params["subset"]["do"]:
        # subset train by threshold
        train = subset_func(train, params["subset"])
        logging.info(f"\nExtracted {train.time.size} images from train.")
        
    train_x = train["anomaly"].values
    theta  = train["params"].values

    logging.info(f"Loaded anomaly training data of shape {train_x.shape}.")
    logging.info(f"Loaded parameters of shape {theta.shape}.")

    # check range of uniform samples
    epsilon = 1e-6
    if not ((compare_u.max() < 1.) and (compare_u.min() > 0.)):
        logging.error("Training uniform samples not in [0,1] range")
        logging.error(f"Max: {compare_u.max()}, Min: {compare_u.min()}")
        compare_u = np.clip(compare_u, epsilon, 1 - epsilon)
    
    if not ((images_u.max() < 1.) and (images_u.min() > 0.)):
        logging.error("Generated uniform samples not in [0,1] range")
        logging.error(f"Max: {images_u.max()}, Min: {images_u.min()}")
        images_u = np.clip(images_u, epsilon, 1 - epsilon)

    # transform images to original scale using invPIT
    distns = [field['distn'] for field in params.fields.values()]
    two_tailed = [field['two_tailed'] for field in params.fields.values()]
    images_x = stats.invPIT(images_u, train_x, theta=theta, distns=distns, two_tailed=two_tailed)
    compare_x = stats.invPIT(compare_u, train_x, theta=theta, distns=distns, two_tailed=two_tailed)

    # save to NetCDF
    data = xr.Dataset(
        {
            "anomaly": (("time", "lat", "lon", "field"), images_x),
            "standardised": (("time", "lat", "lon", "field"), images_y),
            "uniform": (("time", "lat", "lon", "field"), images_u),
            "params": (("lat", "lon", "param", "field"), theta),
        },
        coords={
            "param": [
                "loc_upper", "scale_upper", "shape_upper",
                "loc_lower", "scale_lower", "shape_lower"
                ],
            "field": list(params.fields.keys()),
            "time": (("time"), np.arange(images_x.shape[0])),
            "lat": (("lat"), train["lat"].values),
            "lon": (("lon"), train["lon"].values),
        },
    )

    data.to_netcdf(output.netcdf, format="NETCDF4")
    logging.info(f"Saved generated data to {output.netcdf}.")
    logging.info("Done.")

    # save the comparison dataset
    compare_data = xr.Dataset(
        {
            "anomaly": (("time", "lat", "lon", "field"), compare_x),
            "standardised": (("time", "lat", "lon", "field"), compare_y),
            "uniform": (("time", "lat", "lon", "field"), compare_u),
            "params": (("lat", "lon", "param", "field"), theta),
        },
        coords={
            "param": [
                "loc_upper", "scale_upper", "shape_upper",
                "loc_lower", "scale_lower", "shape_lower"
                ],
            "field": list(params.fields.keys()),
            "time": (("time"), np.arange(compare_y.shape[0])),
            "lat": (("lat"), train["lat"].values),
            "lon": (("lon"), train["lon"].values),
        },
    )

    compare_data.to_netcdf(output.train, format="NETCDF4")
    logging.info(f"Saved training data to {output.train}.")


if __name__     == "__main__":
    # configure logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    main(input, output, params)
    logging.info("Process generated script executed successfully.")