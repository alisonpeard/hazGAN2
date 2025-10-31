"""Load all generated PNGs, inverse-PIT to original scales, and save to NetCDF."""
import os
import logging
import zipfile
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


def load_npy_from_zip(path:str, K=255) -> np.array:
    with zipfile.ZipFile(path, 'r') as zip_ref:
        flist = [f for f in zip_ref.namelist() if f.endswith('.npy')]
        logging.info(f"Found {len(flist)} training images in {input.training_dir}.")
        samples = []
        for f in sorted(flist):
            with zip_ref.open(f) as file:
                sample = np.load(file)
                samples.append(sample / K)
    samples = np.stack(samples, axis=0)
    return samples


def transform_to_scale(y:np.array, image_stats:dict) -> np.array:
    """Rescale samples using image statistics."""
    image_minima = image_stats['min']
    image_maxima = image_stats['max']
    image_range  = image_maxima - image_minima
    logging.info(f"Image statistics: min {image_minima}, max {image_maxima}")
    y_rescaled = y * image_range + image_minima
    y_rescaled = np.flip(y_rescaled, axis=1)
    return y_rescaled


def transform_to_uniform(y:np.array, domain:str) -> np.array:
    inv_transform = getattr(stats, f"inv_{domain}")
    u = inv_transform(y)
    return u


def check_range(u:np.array, epsilon=1e-6):
    if not ((u.max() < 1.) and (u.min() > 0.)):
        logging.error("Generated uniform samples not in [0,1] range")
        logging.error(f"Max: {u.max()}, Min: {u.min()}")
        u = np.clip(u, epsilon, 1 - epsilon)
    return u

def transform_to_physical(u:np.array, ref:xr.Dataset, subset:dict, fields:dict) -> np.array:
    u = check_range(u)
    if subset["do"]:
        ref = subset_func(ref, subset)
        logging.info(f"\nExtracted {ref.time.size} images from train.")
    x_ref = ref["anomaly"].values
    theta  = ref["params"].values
    logging.info(f"Loaded anomaly training data of shape {x_ref.shape}.")
    logging.info(f"Loaded parameters of shape {theta.shape}.")
    distns = [field['distn'] for field in fields.values()]
    two_tailed = [field['two_tailed'] for field in fields.values()]
    x = stats.invPIT(u, x_ref, theta=theta, distns=distns, two_tailed=two_tailed)
    return x


def construct_dataset(images_x, images_y, images_u, theta, lats, lons, fields):
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
            "field": fields,
            "time": (("time"), np.arange(images_x.shape[0])),
            "lat": (("lat"), lats),
            "lon": (("lon"), lons),
        },
    )
    return data


def main(input, output, params):
    # load reference data
    image_stats = np.load(input.image_stats)
    ref = xr.open_dataset(input.training_data)
    theta = ref["params"].values
    lats = ref["lat"].values
    lons = ref["lon"].values

    # load generated images from zipfile
    y = load_npy_from_zip(input.image_dir)
    logging.info(f"Created generated images ndarray of shape {y.shape}.")
    y = transform_to_scale(y, image_stats)
    u = transform_to_uniform(y, params.domain)
    x = transform_to_physical(u, ref, params.subset, params.fields)
    ds = construct_dataset(x, y, u, theta, lats, lons, list(params.fields.keys()))
    ds.to_netcdf(output.netcdf, format="NETCDF4")
    logging.info(f"Saved generated data to {output.netcdf}.")

    # repeat for comparison data
    y = load_npy_from_zip(input.training_dir)
    y = transform_to_scale(y, image_stats)
    u = transform_to_uniform(y, params.domain)
    x = transform_to_physical(u, ref, params.subset, params.fields)
    ds = construct_dataset(x, y, u, theta, lats, lons, list(params.fields.keys()))
    ds.to_netcdf(output.train, format="NETCDF4")
    logging.info(f"Saved training data to {output.train}.")

    logging.info("Done.")


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