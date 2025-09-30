"""Process marginals with fitted tail(s) to create GAN training dataset."""
import os
import sys
import logging
from calendar import month_name

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from typing import Tuple

os.environ["USE_PYGEOS"] = "0"


def load_medians(path:str) -> pd.DataFrame:
    """Load monthly medians and prepare for tensor creation."""
    medians = pd.read_parquet(path)
    medians["lat"] = medians["lat"].astype(float)
    medians["lon"] = medians["lon"].astype(float)
    medians["month_id"] = medians["month"].map(lambda x: list(month_name).index(x))
    medians = medians.sort_values(["month_id", "lat", "lon"]).drop(columns=["month_id"])
    return medians


def load_events(
        events_path:str, metadata_path:str, fields:list
        ) -> Tuple[gpd.GeoDataFrame, pd.Series, float]:
    """Load events data and metadata, return GeoDataFrame and event times."""
    # Load events
    df = pd.read_parquet(events_path)
    df.columns = [col.replace(".", "_") for col in df.columns]
    df[f'day_of_{fields[0]}'] = df.groupby('event')[f'time_{fields[0]}'].rank('dense')
    times_in_df = list(set(pd.to_datetime(df[f'time_{fields[0]}']).values))
    
    # Load metadata
    metadata = pd.read_parquet(metadata_path)
    metadata["time"] = pd.to_datetime(metadata["time"])
    metadata = metadata[metadata["time"].isin(times_in_df)].copy()
    event_times = pd.to_datetime(metadata[['event', 'time']].groupby('event').first()['time'].reset_index(drop=True))
    event_sizes = metadata[["event", "event.size"]].groupby("event").mean()["event.size"]
    yearly_rate = metadata['lambda'].iloc[0]
    
    # Create GeoDataFrame
    df["size"] = df["event"].map(event_sizes)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs="EPSG:4326")
    gdf["lat"] = gdf.geometry.y
    gdf["lon"] = gdf.geometry.x
    gdf = gdf.sort_values(["event", "lat", "lon"]).reset_index(drop=True)
    
    return gdf, event_times, yearly_rate


def validate_data(gdf:gpd.GeoDataFrame, fields:list) -> None:
    """Validate ECDF values are in [0, 1]."""
    for field in fields:
        col = f"ecdf_{field}"
        assert 0 <= gdf[col].min() and gdf[col].max() <= 1, f"ECDF {col} not in [0,1]"


def create_parameters(
        gdf:gpd.GeoDataFrame, h:int, w:int, fields:list
        ) -> np.ndarray:
    """Create hxwx6xk parameter array for distribution parameters."""
    params = np.full((h, w, 6, len(fields)), np.nan, dtype=np.float32)
    
    prefixes = ["thresh_", "scale_", "shape_"]
    suffixes = ["", "_lower", "_upper"]
    param_map = {0: 3, 1: 0, 2: 0}  # suffix index to parameter offset
    
    for k, field in enumerate(fields):
        for i, suffix in enumerate(suffixes):
            for j, prefix in enumerate(prefixes):
                col = f"{prefix}{field}{suffix}"
                if col in gdf.columns:
                    values = gdf[[col, "lat", "lon"]].groupby(["lat", "lon"]).mean()[col].values
                    params[:, :, j + param_map[i], k] = values.reshape([h, w])
    
    return params


def create_tensors(
        gdf:gpd.GeoDataFrame, medians:pd.DataFrame, fields:list
        ) -> dict:
    """Create all training tensors."""
    w, h, n = gdf["lon"].nunique(), gdf["lat"].nunique(), gdf["event"].nunique()
    k = len(fields)
    logging.info(f"Found {n} events.")
    
    # Coordinates
    lats, lons = gdf["lat"].unique(), gdf["lon"].unique()
    
    # Main tensors
    anomalies = gdf[fields].values.reshape([n, w, h, k])
    days_of = gdf[[f"day_of_{fields[0]}"]].values.reshape([n, h, w])
    ecdf = gdf[[f"ecdf_{f}" for f in fields]].values.reshape([n, h, w, k])
    scdf = gdf[[f"scdf_{f}" for f in fields]].values.reshape([n, h, w, k])
    num_months = medians["month"].nunique()
    months = medians["month"].unique().tolist()
    medians = medians.sort_values(["month", "lat", "lon"], ascending=[True, False, True])
    medians = medians[fields].values.reshape([num_months, h, w, k])
    print(f"{medians.shape=}")
    print(f"{medians=}")
    
    # Event data
    return_periods = gdf[["event", "event_rp"]].groupby("event").mean().values.reshape(n)
    sizes = gdf[["event", "size"]].groupby("event").mean().values.reshape(n)
    
    # Parameters
    params = create_parameters(gdf, h, w, fields)
    
    return {
        'coords': (lats, lons),
        'tensors': (
            anomalies, days_of, ecdf, scdf, medians, 
            return_periods, sizes, params
            ),
        'months': months
    }


def main(input, output, params):
    """Main processing function."""
    fields = list(params.fields.keys())
    
    logging.info("Loading data...")
    medians = load_medians(input.medians)
    gdf, times, rate = load_events(input.events, input.metadata, fields)
    
    logging.info("Validating and processing...")
    validate_data(gdf, fields)
    data = create_tensors(gdf, medians, fields)
    
    logging.info("Creating dataset...")
    lats, lons = data['coords']
    anomalies = data["tensors"][0]
    ecdf      = data["tensors"][2]
    scdf      = data["tensors"][3]
    medians   = data["tensors"][4]
    return_periods = data["tensors"][5]
    sizes     = data["tensors"][6]
    params    = data["tensors"][7] # overwriting snakemake params here
    
    ds = xr.Dataset({
        'uniform': (['time', 'lat', 'lon', 'field'], scdf),
        'ecdf': (['time', 'lat', 'lon', 'field'], ecdf),
        'anomaly': (['time', 'lat', 'lon', 'field'], anomalies),
        'medians': (['month', 'lat', 'lon', 'field'], medians),
        'event_rp': (['time'], return_periods),
        'duration': (['time'], sizes),
        'params': (['lat', 'lon', 'param', 'field'], params),
    }, coords={
        'lat': lats, 'lon': lons,
        'time': times,
        'field': fields,
        'param': [
            "loc_upper", "scale_upper", "shape_upper",
            "loc_lower", "scale_lower", "shape_lower"
        ],
        'month': data['months']
    }, attrs={
        'CRS': 'EPSG:4326',
        'yearly_freq': rate,
        'fields_info': str(fields)
    })
    
    logging.info(f"Saving dataset {ds.uniform.shape} to {output.data}")
    ds.to_netcdf(output.data)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(snakemake.log.file), logging.StreamHandler(sys.stdout)]
    )
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    logging.info(f"Input: {input}, Output: {output}, Params: {params}")
    main(input, output, params)