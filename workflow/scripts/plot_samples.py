"""Make plots of data samples.

Note: this is before adding back monthly medians.
"""
import os
import ssl
import numpy as np
import xarray as xr

from src import funcs
from src import plotting
from src import statistics


ssl._create_default_https_context = ssl._create_stdlib_context


def rank_by_intensity_func(ds:xr.Dataset, subset:dict):
    """Rank the dataset using config-defined function and threshold."""
    func = getattr(funcs, subset["func"])
    args = subset["args"]
    thresh = subset["value"]

    for arg in args:
        ds[arg] = ds.sel(field=arg).anomaly

    intensity = func(ds, *args, dim=["lon", "lat"])
    ranking = np.argsort(intensity)[::-1]

    return ranking


def main(input, output, params):
    
    os.makedirs(output.outdir, exist_ok=True)

    transform = getattr(statistics, params.domain)

    # load data and samples
    train = xr.open_dataset(input.train)
    train_x = train["anomaly"].values
    train_u = train["uniform"].values
    train_y = transform(train_u)
    # train_y = train["standardised"].values # if using comparison data instead

    gener = xr.open_dataset(input.generated)
    gener_x = gener["anomaly"].values
    gener_u = gener["uniform"].values
    # gener_y = transform(gener_u) # if using comparison data instead
    gener_y = gener["standardised"].values

    if params.shuffle:
        train_ids = np.random.permutation(train_u.shape[0])
        gener_ids = np.random.permutation(gener_u.shape[0])
    else:
        # sort by intensity fn
        if params["intensity"]["do"]:
            train_ids = rank_by_intensity_func(train, params.intensity)
            gener_ids = rank_by_intensity_func(gener, params.intensity)
        else:
            # just sort by first field
            train_ids = np.argsort(np.max(train_u[..., 0], axis=(1, 2)))[::-1]
            gener_ids = np.argsort(np.max(gener_u[..., 0], axis=(1, 2)))[::-1]

    xmin = params["lon_min"]
    xmax = params["lon_max"]
    ymin = params["lat_min"]
    ymax = params["lat_max"]
    extent = [xmin, xmax, ymin, ymax]

    # make the samples
    for i, field in enumerate(params.fields.keys()):
        metric = params.fields[field].get("units", "")
        cmap   = params.fields[field].get("cmap", "viridis")
        cmap   = plotting.eval_cmap_str(cmap)

        figa = plotting.samples.plot(
            gener_y[gener_ids], train_y[train_ids], field=i, title="",
            extent=extent,
            cbar_label="", cmap=cmap, ndecimals=0
            )
        figb = plotting.samples.plot(
            gener_u[gener_ids], train_u[train_ids], field=i, title="",
            extent=extent,
            cbar_label="", cmap=cmap, ndecimals=1
            )
        figc = plotting.samples.plot(
            gener_x[gener_ids], train_x[train_ids], field=i, title="",
            extent=extent,
            cbar_label=metric, cmap=cmap, alpha=1e-6
            );

        figa.savefig(os.path.join(
            output.outdir, f"{field}_{params.domain}.png"
            ), dpi=300, bbox_inches="tight")
        figb.savefig(os.path.join(
            output.outdir, f"{field}_uniform.png"
            ), dpi=300, bbox_inches="tight")
        figc.savefig(os.path.join(
            output.outdir, f"{field}_anomaly.png"
            ), dpi=300, bbox_inches="tight")

        dates = {
            "generated": gener.time.values[gener_ids],
            "training": train.time.values[train_ids],
        }
        # save to .txt file
        with open(os.path.join(output.outdir, f"{field}_dates.txt"), "w") as f:
            for data in dates["training"]:
                f.write(f"Training date: {data}\n")


if __name__ == "__main__":
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params

    main(input, output, params)