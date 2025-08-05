import os
import numpy as np
import xarray as xr
import logging

from src import funcs
from src import plotting


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

    os.makedirs(output.dir0, exist_ok=True)
    os.makedirs(output.dir1, exist_ok=True)

    # load data and samples
    train = xr.open_dataset(input.train)

    if params["subset"]["do"]:
        # subset train by threshold
        train = subset_func(train, params["subset"])
        logging.info(f"\nExtracted {train.time.size} images from train.")

    gener = xr.open_dataset(input.generated)

    train = train[params.domain].values
    gener = gener[params.domain].values
    
    xmin = params["lon_min"]
    xmax = params["lon_max"]
    ymin = params["lat_min"]
    ymax = params["lat_max"]
    extent = [xmin, xmax, ymin, ymax]

    for _fields in [[0, 1], [0, 2], [1, 2]]:

        logging.info(f"Plotting fields {_fields[0]} and {_fields[1]}")

        fig_smith = plotting.fields.plot(
            gener, train, plotting.fields.smith1990, fields=_fields, figsize=.6, title="",
            cbar_label=r"$\hat\theta$", cmap="Spectral", vmin=1, vmax=4, extent=extent
            )
        fig_smith.savefig(os.path.join(output.dir0, f"smith1990_{_fields[0]}-{_fields[1]}.png"), dpi=300)

        fig_pears = plotting.fields.plot(
            gener, train, plotting.fields.pearson, fields=_fields, figsize=.6,
            title="", cbar_label=r"$r$", vmin=-1, vmax=1, cmap="Spectral_r", extent=extent
            )
        fig_pears.savefig(os.path.join(output.dir0, f"pearson_{_fields[0]}-{_fields[1]}.png"), dpi=300)
    
    # - - - - - Spatial - - - - - - - - 
    logging.info("Resampling to 16 x 16")
    outres = params.outres
    ntrain, inx, _, _ = train.shape
    ngener, _, _, _ = gener.shape
    bin_size = inx // outres

    train_sample = train
    gener_sample = gener

    train_sample = train_sample.reshape(ntrain, outres, bin_size, outres, bin_size, 3).max(3).max(1)
    gener_sample = gener_sample.reshape(ngener, outres, bin_size, outres, bin_size, 3).max(3).max(1)
    logging.info(f"Train shape: {train_sample.shape}")
    logging.info(f"Generated shape: {gener_sample.shape}")

    for _field in [0, 1, 2]:
        logging.info(f"Plotting field {_field}")
        
        # plot Pearson correlation
        fig_pears = plotting.spatial.plot(
            gener_sample, train_sample, plotting.spatial.pearson, field=_field,
            figsize=.6, title="", cbar_label=r"$r$", vmin=-1, vmax=1, cmap="Spectral_r"
            )
        fig_pears.savefig(os.path.join(output.dir1, f"pearson_{_field}.png"), dpi=300)

        logging.info(f"Plotting Smith (1990) field {_field}")
        fig_smith = plotting.spatial.plot(
            gener_sample , train_sample , plotting.spatial.smith1990, field=_field,
            figsize=.6, title="", cbar_label=r"$\hat\theta$", cmap="Spectral",
            vmin=1, vmax=4
            )
        fig_smith.savefig(os.path.join(output.dir1, f"smith1990_{_field}.png"), dpi=300)


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