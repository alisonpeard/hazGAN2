"""Create RGB images from the transformed xarray dataset."""
import os
import numpy as np
import xarray as xr
from PIL import Image
import matplotlib.pyplot as plt
from PIL import Image
import zipfile
import numpy as np
import logging


from src import funcs
import src.python.statistics as stats


def apply_colormap(grayscale_array, colormap_name='Spectral_r'):
    normalized = grayscale_array.astype(float) / 255
    colormap = plt.get_cmap(colormap_name)
    colored = colormap(normalized)
    rgb_uint8 = np.uint8(colored[..., :3] * 255)
    return rgb_uint8


def subset_func(ds:xr.Dataset, subset:dict):
    """Subset the dataset using config-defined function and threshold."""
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
    # load data
    ds = xr.open_dataset(input.data)
    assert ds.uniform.shape[1:] == (params.resx, params.resy, 3), f"Unexpected shape: {ds.uniform.shape}"

    # Make a more extreme dataset
    if params["subset"]["do"]:
        ds = subset_func(ds, params["subset"])
        logging.info(f"\nExtracted {ds.time.size} images.")

    # make PNGs of stacked percentiles
    os.makedirs(output.outdir, exist_ok=True)

    nimgs = ds.time.size
    array = ds.uniform.values
    array = np.flip(array, axis=1) # flip latitude

    if not ((array.max() <= 1.) and (array.min() >= 0.)):
            raise ValueError("Percentiles not in [0,1] range")

    assert array.shape[1:] == (params.resx, params.resy, 3), f"Unexpected shape: {array.shape}"

    # rescale to (0, 1) if domain is not uniform
    if params.domain is not None:
        ppf = getattr(stats, params.domain)
        array = np.clip(array, params.eps, 1-params.eps) # Avoid log(0)
        array = ppf(array)

        # scaling with explicit max / min
        array_min = ppf(1 / params.rp_max)
        array_max = ppf(1 - 1 / params.rp_max)
        logging.info(f"Using {params.domain} with min {array_min} and max {array_max}")

        array = (array - array_min) / (array_max - array_min)

        logging.info("Range: {}--{}".format(array.min(), array.max()))
        logging.info("Shape: {},{}".format(array_min.shape, array_max.shape))

        np.savez(output.image_stats, min=array_min, max=array_max)#, n=n)

    # convert images to RGB and save
    for i in range(nimgs):
        arr = array[i]
        # PNG output
        # arr = np.uint8(arr * 255)
        # img = Image.fromarray(arr, 'RGB')
        # output_path = os.path.join(output.outdir, f"footprint{i}.png")
        # img.save(output_path)

        # npy output
        output_path = os.path.join(output.outdir, f"footprint{i}.npy")
        np.save(output_path, arr * 255)

    # verify saved image
    test_load = np.load(output_path)
    logging.info(f"Saved {nimgs} npy files to {output.outdir}")

    # save to zipfile
    with zipfile.ZipFile(output.zipfile, 'w') as zipf:
        for root, dirs, files in os.walk(output.outdir):
            for file in files:
                if file.endswith('.npy'):
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, os.path.relpath(file_path, output.outdir))
    
    logging.info(f"Saved {output.zipfile} with {len(files)} images")

    

if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
        )
    
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params

    main(input, output, params)