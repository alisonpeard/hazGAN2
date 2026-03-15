"""Create RGB images from the transformed xarray dataset."""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
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
    ds = xr.open_dataset(input.data)
    assert ds.uniform.shape[1:] == (params.resx, params.resy, 3), \
        f"Unexpected shape: {ds.uniform.shape}"

    if params["subset"]["do"]:
        # Make a more extreme dataset
        ds = subset_func(ds, params["subset"])
        logging.info(f"\nExtracted {ds.time.size} events.")

    os.makedirs(output.outdir, exist_ok=True)
    
    nimgs = ds.time.size
    u = ds.uniform.values

    if not ((u.max() < 1.) and (u.min() > 0.)):
        raise ValueError("Percentiles not in (0,1) range")

    # transform to reduced variate
    ppf = getattr(stats, params.domain)
    y = ppf(u)

    # scale back to (0, 1) with rp-based scaling
    ymin = ppf(1 / params.rpmax)
    ymax = ppf(1 - 1 / params.rpmax)

    logging.info(f"Using {params.domain} with min {ymin} and max {ymax}")
    assert ymin < y.min() < y.max() < ymax, \
        f"Data outside expected range. Check rpmax/domain.\n" \
        f"Data range: {y.min():.4f} - {y.max():.4f}.\n" \
        f"Expected range: {ymin:.4f} - {ymax:.4f}.\n" \
        f"Current rpmax: {params.rpmax:,f}.\n" \
        f"Uniform range: {u.min():.4f} - {u.max():.4f}. " 
    
    y_scaled = (y - ymin) / (ymax - ymin)
    logging.info("Range: {} - {}".format(y_scaled.min(), y_scaled.max()))
    
    # save image stats for inverse scaling
    np.savez(output.image_stats, min=ymin, max=ymax)
    logging.info(f"Saved image stats to {output.image_stats}")

    # convert images to RGB and save
    for i in range(nimgs):
        output_path = os.path.join(output.outdir, f"footprint{i}.npy")
        np.save(output_path, y_scaled[i] * 255)

    # test load saved image
    _ = np.load(output_path)
    logging.info(f"Saved {nimgs} npy files to {output.outdir}")

    # zip everything for easier file transfer
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