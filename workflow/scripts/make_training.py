"""
Process the marginals with fitted tails to make training dataset for the GAN.

NOTE: generic names (shape, scale, loc) will be used regardless of distribution.
These should be set to NaN accordingly.

TODO: This script can be made way more concise.
"""
# %%
import os
os.environ["USE_PYGEOS"] = "0"
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import logging

from calendar import month_name as month


if __name__ == "__main__":
    # configure logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    
    # load parameters
    DATA     = snakemake.input.data_all
    EVENTS   = snakemake.input.events
    METADATA = snakemake.input.metadata
    # MEDIANS  = snakemake.input.medians # TODO
    OUTPUT   = snakemake.output.data
    FIELDS = [key for key in snakemake.params.fields.keys()]

    # load coordinates
    data = xr.open_dataset(DATA)
    coords = data['grid'].to_dataframe().reset_index()
    coords = gpd.GeoDataFrame(
        coords, geometry=gpd.points_from_xy(coords['lon'], coords['lat'])
        ).set_crs("EPSG:4326")
    monthly_medians = data.groupby("time.month").median(dim="time").to_dataframe().reset_index()
    monthly_medians = monthly_medians[['grid', 'month'] + FIELDS]

    # load fitted data                                                                                                        
    df = pd.read_parquet(EVENTS)
    df = df.merge(coords, on="grid")
    df.columns = [col.replace(".", "_") for col in df.columns]
    df[f'day_of_{FIELDS[0]}'] = df.groupby('event')[f'time_{FIELDS[0]}'].rank('dense')

    # get event times and durations
    events = pd.read_parquet(METADATA)
    times = pd.to_datetime(
        events[['event', 'time']]
        .groupby('event').first()['time']
        .reset_index(drop=True)
        )
    rate = events['lambda'][0]
    events = events[["event", "event.size", "lambda"]].groupby("event").mean()
    events = events.to_dict()["event.size"]
    df["size"] = df["event"].map(events)
    gdf = gpd.GeoDataFrame(df, geometry="geometry").set_crs("EPSG:4326")

    # check ecdfs are in (0, 1)
    ecdf_cols = [f"ecdf_{field}" for field in FIELDS]
    for col in ecdf_cols:
        assert gdf[col].max() <= 1, f"ECDF values for {col} should be <= 1"
        assert gdf[col].min() >= 0, f"ECDF values for {col} should be >= 0"

    # merge in monthly medians
    logging.warning("Calculating medians on the fly, this should be done earlier.")
    nduplicates = monthly_medians.groupby(['month', 'grid']).count().max().max()
    assert nduplicates == 1, "Monthly medians not unique"
    del nduplicates
    monthly_medians = monthly_medians.groupby(["month", "grid"] + FIELDS).mean().reset_index()
    monthly_medians["month"] = monthly_medians["month"].apply(lambda x: month[x])

    for field in FIELDS:
        gdf[f"month_{field}"] = pd.to_datetime(gdf[f"time_{field}"]).dt.month.map(lambda x: month[x])
        n = len(gdf)
        gdf = gdf.join(
            monthly_medians[['month', 'grid', field]].set_index(["month", "grid"]),
            on=[f"month_{field}", "grid"],
            rsuffix="_median"
            )
        assert n == len(gdf), "Merge failed"
        del gdf[f'month_{field}']

    # use lat and lon columns to label grid points in (i,j) format
    gdf["lat"] = gdf["geometry"].apply(lambda x: x.y)
    gdf["lon"] = gdf["geometry"].apply(lambda x: x.x)
    gdf = gdf.sort_values(["lat", "lon", "event"], ascending=[True, True, True])

    #  make netcdf file
    nfields = len(FIELDS)
    nx      = gdf["lon"].nunique()
    ny      = gdf["lat"].nunique()
    T       = gdf["event"].nunique()

    # make training tensors
    gdf  = gdf.sort_values(["event", "lat", "lon"], ascending=[True, True, True]) # [T, i, j, field]
    grid = gdf["grid"].unique().reshape([ny, nx])
    lat  = gdf["lat"].unique()
    lon  = gdf["lon"].unique()
    X    = gdf[FIELDS].values.reshape([T, ny, nx, nfields])
    D    = gdf[[f"day_of_{FIELDS[0]}"]].values.reshape([T, ny, nx])
    U0   = gdf[[f"ecdf_{field}" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    U1   = gdf[[f"scdf_{field}" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    M    = gdf[[f"{field}_median" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    z    = gdf[["event", "event_rp"]].groupby("event").mean().values.reshape(T)
    s    = gdf[["event", "size"]].groupby("event").mean().values.reshape(T)

    # TODO: add checks here of lifetime max / total
    logging.warning("Lifetime max / total not checked")
    
    gpd_params = ([f"thresh_{var}" for var in FIELDS] + [f"scale_{var}" for var in FIELDS] + [f"shape_{var}" for var in FIELDS])
    gdf_params = (gdf[[*gpd_params, "lon", "lat"]].groupby(["lat", "lon"]).mean().reset_index())
    thresh = np.array(gdf_params[[f"thresh_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
    scale = np.array(gdf_params[[f"scale_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
    shape = np.array(gdf_params[[f"shape_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
    params = np.stack([thresh, scale, shape], axis=-2)
    logging.info("Parameters shape: {}".format(params.shape))

    # make an xarray dataset for training
    ds = xr.Dataset({'uniform': (['time', 'lat', 'lon', 'field'], U1),
                    'ecdf': (['time', 'lat', 'lon', 'field'], U0),
                    'anomaly': (['time', 'lat', 'lon', 'field'], X),
                    'medians': (['time', 'lat', 'lon', 'field'], M),
                    f'day_of_{FIELDS[0]}': (['time', 'lat', 'lon'], D),
                    'event_rp': (['time'], z),
                    'duration': (['time'], s),
                    'params': (['lat', 'lon', 'param', 'field'], params),
                    'grid': (['lat', 'lon'], grid),
                    },
                    coords={'lat': (['lat'], lat),
                            'lon': (['lon'], lon),
                            'time': times,
                            'field': FIELDS,
                            'param': ['loc', 'scale', 'shape']
                    },
                    attrs={'CRS': 'EPSG:4326',
                            'yearly_freq': rate,
                            'fields_info': str(FIELDS)})

    # save
    logging.info("Finished! Saving to netcdf...")
    ds.to_netcdf(OUTPUT)
    logging.info("Saved to{}".format(OUTPUT))

    