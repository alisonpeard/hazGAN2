# %%
import os
import numpy as np
import xarray as xr
from PIL import Image
import matplotlib.pyplot as plt
from PIL import Image
import zipfile
import numpy as np
import logging


def apply_colormap(grayscale_array, colormap_name='Spectral_r'):
    normalized = grayscale_array.astype(float) / 255
    colormap = plt.get_cmap(colormap_name)
    colored = colormap(normalized)
    rgb_uint8 = np.uint8(colored[..., :3] * 255)
    return rgb_uint8

if __name__ == "__main__":
    # TODO: 1. make extreme extraction optional
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
        )
    DATA     = snakemake.input["data"]
    OUTDIR   = snakemake.output["outdir"]
    ZIPFILE  = snakemake.output["zipfile"]
    STATS    = snakemake.output["image_stats"]
    THRESH   = snakemake.params["event_subset"]
    DO_SUBSET = snakemake.params["do_subset"]
    DOMAIN   = snakemake.params["domain"]
    EPS      = snakemake.params["eps"]
    RESX     = snakemake.params["resx"]
    RESY     = snakemake.params["resy"]

    # load data
    ds = xr.open_dataset(DATA)

    # Make a more extreme dataset
    if DO_SUBSET:
        ds['intensity'] = ds.sel(field=THRESH["field"]).anomaly.apply(THRESH["func"], dims=['lon', 'lat'])
        mask = (ds['intensity'] > THRESH["value"]).values
        idx  = np.where(mask)[0]
        ds   = ds.isel(time=idx)

        logging.info(f"\nFound {ds.time.size} images with {THRESH['func']} {THRESH['field']} exceeding {THRESH['value']} m/s")

    # make PNGs of stacked percentiles
    os.makedirs(OUTDIR, exist_ok=True)

    nimgs = ds.time.size
    array = ds.uniform.values
    array = np.flip(array, axis=1) # flip latitude

    if not ((array.max() <= 1.) and (array.min() >= 0.)):
            raise ValueError("Percentiles not in [0,1] range")

    assert array.shape[1:] == (RESX, RESY, 3), f"Unexpected shape: {array.shape}"

    # rescale to (0, 1) if data is Gumbel distributed
    if DOMAIN == "gumbel":
        array = np.clip(array, EPS, 1-EPS) # Avoid log(0)
        array = -np.log(-np.log(array))
        array_min = np.min(array, axis=(0, 1, 2), keepdims=True)
        array_max = np.max(array, axis=(0, 1, 2), keepdims=True)
        n = len(array)

        # scale to (0, 1)
        array = (array - array_min) / (array_max - array_min)
        array = (array * (n - 1) + 1) / (n + 1)
        # NOTE: original = ((scaled * (n+1) - 1) / (n-1)) * (max - min) + min

        logging.info("Range: {}--{}".format(array.min(), array.max()))
        logging.info("Shape: {},{}".format(array_min.shape, array_max.shape))
        np.savez(STATS, min=array_min, max=array_max, n=n)

    # convert images to RGB and save
    for i in range(nimgs):
        arr = array[i]
        arr = np.uint8(arr * 255)
        img = Image.fromarray(arr, 'RGB')
        output_path = os.path.join(OUTDIR, f"footprint{i}.png")
        img.save(output_path)

    # verify saved image
    test_load = Image.open(output_path)
    logging.info(f"Saved {nimgs} PNG images to {OUTDIR}")

    # save to zipfile
    with zipfile.ZipFile(ZIPFILE, 'w') as zipf:
        for root, dirs, files in os.walk(OUTDIR):
            for file in files:
                if file.endswith('.png'):
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, os.path.relpath(file_path, OUTDIR))
    
    logging.info(f"Saved {ZIPFILE} with {len(files)} images")

    