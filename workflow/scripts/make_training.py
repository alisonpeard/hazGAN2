"""
Process the marginals with fitted tails to make training dataset for the GAN.
"""
# %%
import os
os.environ["USE_PYGEOS"] = "0"
import pandas as pd
import geopandas as gpd
import xarray as xr
import logging

from calendar import month_name as month

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

FIELDS = [key for key in snakemake.params.fields.keys()]

# for snakemake in future
INFILES  = snakemake.input
OUTFILES = snakemake.output

if __name__ == "__main__":
    # load coordinates
    coords = xr.open_dataset(INFILES[0])
    coords = coords['grid'].to_dataframe().reset_index()
    coords = gpd.GeoDataFrame(
        coords, geometry=gpd.points_from_xy(coords['lon'], coords['lat'])
        ).set_crs("EPSG:4326")

    # load fitted data                                                                                                        
    df = pd.read_parquet(INFILES[1])

    df = df.merge(coords, on="grid")
    df.columns = [col.replace(".", "_") for col in df.columns]
    df['day_of_event'] = df.groupby('event')['time_u10'].rank('dense')

    #  add event durations
    events = pd.read_parquet(INFILES[2])
    rate = events['lambda'][0]
    events = events[["event", "event.size", "lambda"]].groupby("event").mean()
    events = events.to_dict()["event.size"]
    df["size"] = df["event"].map(events)
    gdf = gpd.GeoDataFrame(df, geometry="geometry").set_crs("EPSG:4326")

    # load event time data
    events = pd.read_parquet(INFILES[2])
    times = pd.to_datetime(events[['event', 'time']].groupby('event').first()['time'].reset_index(drop=True))
    
    
    # return gdf # TODO: remove this line later
    #  important: check ecdfs are in (0, 1)
    ecdf_cols = [f"ecdf_{field}" for field in FIELDS]
    for col in ecdf_cols:
        assert gdf[col].max() <= 1, f"ECDF values for {col} should be <= 1"
        assert gdf[col].min() >= 0, f"ECDF values for {col} should be >= 0"

    # merge in monthly medians
    monthly_medians = pd.read_csv(INFILES[3], index_col="month")
    assert monthly_medians.groupby(['month', 'grid']).count().max().max() == 1, "Monthly medians not unique"
    monthly_medians = monthly_medians.groupby(["month", "grid"]).mean().reset_index()

    for field in FIELDS:
        gdf[f"month_{field}"] = pd.to_datetime(gdf[f"time_{field}"]).dt.month.map(lambda x: month[x])
        n = len(gdf)
        gdf = gdf.join(monthly_medians[['month', 'grid', field]].set_index(["month", "grid"]), on=[f"month_{field}", "grid"], rsuffix="_median")
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
    D    = gdf[["day_of_event"]].values.reshape([T, ny, nx])
    U0   = gdf[[f"ecdf_{field}" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    U1   = gdf[[f"scdf_{field}" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    M    = gdf[[f"{field}_median" for field in FIELDS]].values.reshape([T, ny, nx, nfields])
    z    = gdf[["event", "event_rp"]].groupby("event").mean().values.reshape(T)
    s    = gdf[["event", "size"]].groupby("event").mean().values.reshape(T)

    # TODO: add checks here of lifetime max / total

    # parameters for parametric tail fits
    # !TODO resolve how to handle 2/3 parameters with different names
    params = [f"params_{field}" for field in FIELDS]
    params = (gdf[[*params, "lon", "lat"]].groupby(["lat", "lon"]).mean().reset_index())
    # thresh = np.array(gdf_params[[f"thresh_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
    # scale = np.array(gdf_params[[f"scale_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
    # shape = np.array(gdf_params[[f"shape_{var}" for var in FIELDS]].values.reshape([ny, nx, nfields]))
    # params = np.stack([thresh, scale, shape], axis=-2)

    # make an xarray dataset for training
    ds = xr.Dataset({'uniform': (['time', 'lat', 'lon', 'field'], U1),
                    'ecdf': (['time', 'lat', 'lon', 'field'], U0),
                    'anomaly': (['time', 'lat', 'lon', 'field'], X),
                    'medians': (['time', 'lat', 'lon', 'field'], M),
                    'day_of_event': (['time', 'lat', 'lon'], D),
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
                            'u10': '10m Wind Speed [ms-1]',
                            'tp': 'Total Precipitation [m]',
                            'mslp': 'Mean Sea Level Pressure [Pa]',
                            'yearly_freq': rate})

    # save
    logging.info("Finished! Saving to netcdf...")
    ds.to_netcdf(OUTFILES[0])
    logging.info("Saved to", OUTFILES[0])

    