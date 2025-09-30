import logging
import pandas as pd
import numpy as np
import xarray as xr

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

    _dataset = params.dataset.upper()

    # load data and samples
    logging.info(f"Plotting barcharts for {_dataset} dataset, month {params.month}")

    train = xr.open_dataset(input.train)
    medians = train.sel(month=[params.month]).medians.values
    
    if params["event_subset"]["do"]:
        # subset train by threshold
        train = subset_func(train, params["event_subset"])
        logging.info(f"\nExtracted {train.time.size} images from train.")
    
    train_x = train["anomaly"].values

    gener = xr.open_dataset(input.generated)
    gener_x = gener["anomaly"].values

    # add medians for month back to samples
    logging.info(f"{medians.shape=}, {train_x.shape=}, {gener_x.shape=}")
    train_x += medians
    gener_x += medians

    bar_width = 0.25
    fig, ax = plotting.misc.saffirsimpson_barchart(gener_x, train_x, bar_width=bar_width, title="")

    if False: # load and add IBTrACS
        raise NotImplementedError("IBTrACS not implemented")
        from hazGAN import saffirsimpson

        ibtracs = pd.read_csv("/Users/alison/Documents/DPhil/paper1.nosync/training/18x22/ibtracs_dates.csv")
        ibtracs["time"] = pd.to_datetime(ibtracs["time"])
        ibtracs = ibtracs.sort_values("time", ascending=False)
        ibtracs = ibtracs[ibtracs['event'] != "Not_Named"]
        ibtracs_storms = ibtracs.groupby("event").agg({"time": ["min", "max"], "wind": "max"})
        ibtracs_storms.columns = ["start", "end", "wind"]
        stormibtracs_stormss = ibtracs_storms.sort_values("start", ascending=False)
        ibtracs_storms = ibtracs_storms.sort_values("wind", ascending=False)
        nstorms = len(ibtracs_storms)
        ibtracs_storms['Category'] = ibtracs_storms['wind'].apply(saffirsimpson)
        cat_counts = ibtracs_storms['Category'].value_counts().sort_index()
        cat_density = cat_counts / len(ibtracs_storms)
        all_categories = np.arange(-1, 6)
        ibtracs_density = pd.Series(cat_density.get(cat, 0) for cat in all_categories)
        r = np.arange(len(all_categories))
        r = [x + 2 * bar_width + 0.01 for x in r]

        ax.bar(r, ibtracs_density, width=bar_width, color="#F6F5EE", label="IBTrACS",
            edgecolor='#666666', linewidth=0.5)
        
    ax.legend()
    fig.tight_layout()
    fig.savefig(output.figure, dpi=300)


if __name__ == "__main__":
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