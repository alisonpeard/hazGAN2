"""Examine the output of extract_events.R"""
# %%
import os
import numpy as np
import pandas as pd
import xarray as xr
from calendar import month_name
import matplotlib.pyplot as plt

wd = os.path.join("..", "..", "results", "processing")
resources_dir = os.path.join("..", "..", "resources")
data_all = "data_all.nc"
medians  = "medians.parquet"
metadata = "event_metadata.parquet"
daily    = "daily.parquet"
storms   = "cyclones_midlands.csv"
fitted   = "events.parquet"


if __name__ == "__main__":
    # load the data
    data_all = xr.open_dataset(os.path.join(wd, data_all))
    storms = pd.read_csv(os.path.join(resources_dir, storms))
    metadata = pd.read_parquet(os.path.join(wd, metadata))
    medians_out = pd.read_parquet(os.path.join(wd, medians))
    daily_out = pd.read_parquet(os.path.join(wd, daily))

    # %% 
    print(metadata["event"].nunique()) # 1592
    print(metadata.shape)  # 12410880 rows × 6 columns

    # %%
    fitted = pd.read_parquet(os.path.join(wd, fitted))
    print(fitted.shape)  # 6086656 rows × 55 columns
    print(fitted["event"].nunique()) # 1485

    # %%


    # %% examine the input storms for NAs
    storms.head()
    """funcs.R::identify_events()

    event_data <- read.csv(input_path)
    event_data <- event_data[event_data$wind > 0.0, ]
    event_data$time <- as.Date(event_data$date)
    event_data$year <- format(event_data$time, "%Y")
    event_data$event <- as.numeric(factor(event_data$cyclone_id))

    """
    storms = storms[storms["wind"] > 0.0].copy()
    storms["time"] = pd.to_datetime(storms["date"])
    storms["year"] = storms["time"].dt.year
    storms["event"] = pd.factorize(storms["cyclone_id"])[0] + 1  # start from 1
    storms.head()

    storms["event"].isna().any() # np.False_

    # %% examine event_metadata.parquet for NAs
    metadata[metadata["event"].isna()]  # 8536064 rows × 6 columns
    metadata[metadata["event"].notna()] # 12410880 rows × 6 columns

    # %% check 1: validate deseasonalisation (will need to redo once I change this)
    columns = ["lat", "lon", "time", "u10", "v10", "r30"]
    sortby = ["lat", "lon", "time"]

    medians_out["lat"] = medians_out["lat"].astype(float)
    medians_out["lon"] = medians_out["lon"].astype(float)
    medians_out["u10"] = medians_out["u10"].astype(np.float32)
    medians_out["v10"] = medians_out["v10"].astype(np.float32)
    medians_out["r30"] = medians_out["r30"].astype(np.float32)
    medians_out = medians_out.sort_values(by=["lat", "lon", "month"]).reset_index(drop=True)

    daily_out["lon"]  = daily_out["lon"].astype(float)
    daily_out["lat"]  = daily_out["lat"].astype(float)
    daily_out["time"] = daily_out["time"].astype("datetime64[ns]")
    daily_out["r30"]  = daily_out["r30"].astype(np.float32)
    daily_out = daily_out[columns].sort_values(by=sortby).reset_index(drop=True)

    if False: # actual test
        daily = data_all.to_dataframe().reset_index()

        daily["time"] = pd.to_datetime(daily["time"])
        daily["month"] = daily["time"].dt.month

        medians = daily[["lat", "lon", "month", "u10", "v10", "r30"]]
        medians = medians.groupby(["lat", "lon", "month"])

        medians_test = medians.median().reset_index()
        medians_test["month"] = medians_test["month"].apply(lambda x: month_name[x])

        medians_test["u10"] = medians_test["u10"].astype(np.float32)
        medians_test["v10"] = medians_test["v10"].astype(np.float32)
        medians_test["r30"] = medians_test["r30"].astype(np.float32)

        medians_test = medians_test.sort_values(by=["lat", "lon", "month"]).reset_index(drop=True)

        pd.testing.assert_frame_equal(
            medians_test, medians_out, atol=1e-6
        )
        del medians_test

        medians = medians.transform("median")
        
        daily["u10"]  -= medians["u10"]
        daily["v10"]  -= medians["v10"]
        daily["r30"] -= medians["r30"]

        del daily["month"]

        # make sure same format as referencedaily = daily[columns].sort_values(by=sortby).reset_index(drop=True)

        pd.testing.assert_frame_equal(
            daily, daily_out, atol=1e-6
            )


    # %% visualise the distributions of the anomalies
    hist_kws = {"color": "lightgrey", "linewidth": 0.5, "edgecolor": "k",
                "bins": 50, "density": True}
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 8))

    daily_out.u10.plot.hist(**hist_kws, xlabel="seasonal anomaly u10 (m/s)", ylabel="Density", ax=axs[0])
    daily_out.v10.plot.hist(**hist_kws, xlabel="seasonal anomaly v10 (m/s)", ylabel="Density", ax=axs[1])
    daily_out.r30.plot.hist(**hist_kws, xlabel="seasonal anomaly r30 (m)", ylabel="Density", ax=axs[2])
    
    plt.tight_layout()

    # %% figure out why some indices are missing
    # Unexpected shape: gdf.shape=(6086656, 55) nfields=3, nx=64, ny=64, T=1485
    # expected number of rows: 6086656
    events_path = os.path.join(wd, "events.parquet")
    events = pd.read_parquet(events_path)

    nevents = events["event"].nunique()
    nx = events["lon"].nunique()
    ny = events["lat"].nunique()
    nfields = 3

    # find out what combinations are missing
    fields = ["u10", "v10", "r30"]
    grouped = events.groupby(["lat", "lon", "event"])
    grouped_first = grouped.first()

    # %% 
    unique_events = events["event"].unique()
    # %%
    # lon: -10.875 -> 5.125
    # lat: 49 -> 64.75
    # events: 1 -> 1485

    import numpy as np
    from itertools import product

    lons = np.linspace(-10.875, 5.125, 64)
    lats = np.linspace(49, 64.75, 64)
    event_list = np.arange(0, 1486)
    fields = ["u10", "v10", "r30"]

    all_combinations = pd.DataFrame(
        list(product(
            lons,
            lats,
            event_list
        )),
        columns=["lat", "lon", "event"]
    ) # ==> 6082560 / 6086656


    # %%
    print(f"nevents: {nevents}, nfields: {nfields}, nx: {nx}, ny: {ny}")
    print(f"Expected len: ({nevents * ny * nx * nfields})")
    print(f"Actual len: {events.shape[0]}")

    # %% 
    events_meta = pd.read_parquet(os.path.join(wd, "event_metadata.parquet"))
    print(events_meta.shape)

    # %%
    daily_out = pd.read_parquet(os.path.join(wd, daily))
    daily_out.head()

    # %%
    input_data = pd.read_parquet(os.path.join(wd, "../../resources/cyclones_midlands.csv"))
    input_data.head()

     # %%
    # check 2: validate storm matching

    # check 3: how big is the difference between storms:wind vs data_all:u10?

    # check 4: distributions of u10 across all storms

    # check 5: distribution of v10 across all storms

    # check 6: distributions of r30 across all storms
    
# %%