"""
Process the marginals with fitted tails to make training dataset for the GAN.

NOTE: generic names (shape, scale, loc) will be used regardless of distribution.
These should be set to NaN accordingly.

TODO: This script is bloated, it can be made way more concise.
"""
# %%
import os
os.environ["USE_PYGEOS"] = "0"
import numpy as np
import pandas as pd
# from shapely.geometry import Point
import geopandas as gpd
import xarray as xr
import logging

from calendar import month_name

if __name__ == "__main__":
    # configure logging
    logging.basicConfig(
        filename=snakemake.log.file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    
    # load parameters
    EVENTS   = snakemake.input.events
    METADATA = snakemake.input.metadata
    MEDIANS  = snakemake.input.medians
    OUTPUT   = snakemake.output.data
    FIELDS = [key for key in snakemake.params.fields.keys()]

    # load monthly medians
    medians = pd.read_parquet(MEDIANS)
    medians["lat"] = medians["lat"].astype(float)
    medians["lon"] = medians["lon"].astype(float)
    msg = f"Shape of monthly medians: {medians.shape=}"
    logging.info(msg) # 49162 is 64 x 64 x 12

    # load fitted data                                                                                                        
    df = pd.read_parquet(EVENTS)
    df.columns = [col.replace(".", "_") for col in df.columns]
    df[f'day_of_{FIELDS[0]}'] = df.groupby('event')[f'time_{FIELDS[0]}'].rank('dense')
    logging.info(f"Latitudes: {df.lat.nunique()=}")

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

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(
        df["lon"],
        df["lat"],
        crs="EPSG:4326"
    )).set_crs("EPSG:4326")

    # check ecdfs are in (0, 1)
    ecdf_cols = [f"ecdf_{field}" for field in FIELDS]
    for col in ecdf_cols:
        assert gdf[col].max() <= 1, f"ECDF values for {col} should be <= 1"
        assert gdf[col].min() >= 0, f"ECDF values for {col} should be >= 0"

    # use lat and lon columns to label grid points in (i,j) format
    gdf["lat"] = gdf["geometry"].apply(lambda x: x.y)
    gdf["lon"] = gdf["geometry"].apply(lambda x: x.x)
    gdf = gdf.sort_values(["lat", "lon", "event"], ascending=[True, True, True])

    #  make netcdf file
    nfields = len(FIELDS)
    nx      = gdf["lon"].nunique()
    ny      = gdf["lat"].nunique()
    T       = gdf["event"].nunique()

    # process monthly medians
    medians["month_id"] = medians["month"].map(lambda x: list(month_name).index(x))
    medians = medians.sort_values(["month_id", "lat", "lon"], ascending=[True, True, True])
    medians = medians.drop(columns=["month_id"])

    # make training tensors
    gdf  = gdf.sort_values(["event", "lat", "lon"], ascending=[True, True, True]) # [T, i, j, field]
    lat  = gdf["lat"].unique()
    lon  = gdf["lon"].unique()
    X    = gdf[FIELDS].values.reshape([T, ny, nx, nfields])
    D    = gdf[[f"day_of_{FIELDS[0]}"]].values.reshape([T, ny, nx])
    U0   = gdf[[f"ecdf_{field}" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    U1   = gdf[[f"scdf_{field}" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    z    = gdf[["event", "event_rp"]].groupby("event").mean().values.reshape(T)
    s    = gdf[["event", "size"]].groupby("event").mean().values.reshape(T)
    M    = medians[FIELDS].values.reshape([12, ny, nx, nfields])

    assert lat.shape == (ny,), f"Unexpected shape: {lat.shape=}"
    assert lon.shape == (nx,), f"Unexpected shape: {lon.shape=}"
    assert X.shape == (T, ny, nx, nfields), f"Unexpected shape: {X.shape=}"
    assert D.shape == (T, ny, nx), f"Unexpected shape: {D.shape=}"
    assert U0.shape == (T, ny, nx, nfields), f"Unexpected shape: {U0.shape=}"
    assert U1.shape == (T, ny, nx, nfields), f"Unexpected shape: {U1.shape=}"

    print(f"Shape {nx=}")
    print(f"Shape {ny=}")
    print(f"Shape {T=}")
    print(f"Shape {nfields=}")
    print(f"Shape {lat.shape=}")
    print(f"Shape {lon.shape=}")
    print(f"Shape {U0.shape=}")
    print(f"Shape {U1.shape=}")

    # TODO: add checks here of lifetime max / total
    logging.warning("Lifetime max / total not checked")
    
    if False: # old
        gpd_params = ([f"thresh_{var}" for var in FIELDS] + [f"scale_{var}" for var in FIELDS] + [f"shape_{var}" for var in FIELDS])
        gdf_params = (gdf[[*gpd_params, "lon", "lat"]].groupby(["lat", "lon"]).mean().reset_index())
        thresh = np.array(gdf_params[[f"thresh_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
        scale = np.array(gdf_params[[f"scale_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
        shape = np.array(gdf_params[[f"shape_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
        params = np.stack([thresh, scale, shape], axis=-2)
        logging.info("Parameters shape: {}".format(params.shape))


    # %% dev
    import pandas as pd
    import numpy as np

    param_cols = [
        "thresh_u10_lower", "thresh_u10_upper", "thresh_r30", "thresh_v10",
        "scale_u10_lower", "scale_u10_upper", "scale_r30", "scale_v10",
        "shape_u10_lower", "shape_u10_upper", "shape_r30", "shape_v10"
        ]
    
    FIELDS = ["u10", "v10", "r30"]

    param_prefixes = ["thresh_", "scale_", "shape_"]
    param_suffixes = ["", "_lower", "_upper"]


    # make params H x W x 6 x K of np.nan
    ny, nx = 64, 64
    nfields = len(FIELDS)
    params = np.full((ny, nx, 6, nfields), np.nan, dtype=np.float32)


    gdf = np.random.random((ny*nx, len(param_cols) + 2))
    gdf = pd.DataFrame(gdf, columns=["lat", "lon"] + param_cols)

    for k, field in enumerate(FIELDS):
        for i, suffix in enumerate(param_suffixes):
            for j, prefix in enumerate(param_prefixes):
                param_col = f"{prefix}{field}{suffix}"
                print(f"{i=}, {j=}, {param_col=}, {param_col in gdf.columns=}")
                if param_col in gdf.columns:

                    gdf_param = gdf[[param_col] + ["lat", "lon"]]
                    gdf_param = gdf_param.groupby(["lat", "lon"]).mean().reset_index()
                    param = gdf_param[param_col].values.reshape([ny, nx])

                    
                    _i = 3 * [0, 1, 0][i]
                    print(f"{_i=}, {_i+j=}, {params.shape=}")
                    params[:, :, j + _i, k] = param

    params


    # %%

    # make an xarray dataset for training
    ds = xr.Dataset({'uniform': (['time', 'lat', 'lon', 'field'], U1),
                    'ecdf': (['time', 'lat', 'lon', 'field'], U0),
                    'anomaly': (['time', 'lat', 'lon', 'field'], X),
                    'medians': (['month', 'lat', 'lon', 'field'], M),
                    f'day_of_{FIELDS[0]}': (['time', 'lat', 'lon'], D),
                    'event_rp': (['time'], z),
                    'duration': (['time'], s),
                    'params': (['lat', 'lon', 'param', 'field'], params),
                    },
                    coords={'lat': (['lat'], lat),
                            'lon': (['lon'], lon),
                            'time': times,
                            'field': FIELDS,
                            'param': ['loc', 'scale', 'shape'],
                            'month': medians["month"].unique()
                    },
                    attrs={'CRS': 'EPSG:4326',
                            'yearly_freq': rate,
                            'fields_info': str(FIELDS)})

    assert ds.uniform.shape[1:] == (ny, nx, nfields), f"Unexpected shape: {ds.uniform.shape=}"
    logging.info(f"Shape of dataset: {ds.uniform.shape=}")
    print(f"Shape of dataset: {ds.uniform.shape=}")
    # save
    logging.info("Finished! Saving to netcdf...")
    ds.to_netcdf(OUTPUT)
    logging.info("Saved to {}".format(OUTPUT))

    
# %%
# %%