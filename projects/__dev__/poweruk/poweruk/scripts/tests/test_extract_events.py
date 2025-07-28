"""Examine the output of extract_events.R"""
# %%
import os
import numpy as np
import pandas as pd
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

wd = os.path.join("..", "..", "results", "processing")
resources_dir = os.path.join("..", "..", "resources")
data_all = "data_all.nc"
medians  = "medians.parquet"
metadata = "event_metadata.parquet"
daily    = "daily.parquet"
storms   = "cyclones_midlands.csv"


if __name__ == "__main__":
    # load the data
    data_all = xr.open_dataset(os.path.join(wd, data_all))
    storms = pd.read_csv(os.path.join(resources_dir, storms))
    metadata = pd.read_parquet(os.path.join(wd, metadata))
    medians_out = pd.read_parquet(os.path.join(wd, medians))
    daily_out = pd.read_parquet(os.path.join(wd, daily))

    # %% test 1: validate deseasonalisation
    daily = data_all.to_dataframe().reset_index()

    daily["time"] = pd.to_datetime(daily["time"])
    daily["month"] = daily["time"].dt.month

    medians = daily[["lat", "lon", "month", "vx", "dx", "r30"]]
    medians = medians.groupby(["lat", "lon", "month"])

    from calendar import month_name

    medians_test = medians.median().reset_index()
    medians_test["month"] = medians_test["month"].apply(lambda x: month_name[x])

    medians_out["lat"] = medians_out["lat"].astype(float)
    medians_out["lon"] = medians_out["lon"].astype(float)
    medians_out["vx"] = medians_out["vx"].astype(np.float32)
    medians_out["dx"] = medians_out["dx"].astype(np.float32)
    medians_out["r30"] = medians_out["r30"].astype(np.float32)

    medians_test["vx"] = medians_test["vx"].astype(np.float32)
    medians_test["dx"] = medians_test["dx"].astype(np.float32)
    medians_test["r30"] = medians_test["r30"].astype(np.float32)

    medians_test = medians_test.sort_values(by=["lat", "lon", "month"]).reset_index(drop=True)
    medians_out = medians_out.sort_values(by=["lat", "lon", "month"]).reset_index(drop=True)

    pd.testing.assert_frame_equal(
        medians_test, medians_out, atol=1e-6
    )
    del medians_test

    medians = medians.transform("median")
    
    daily["vx"]  -= medians["vx"]
    daily["dx"]  -= medians["dx"]
    daily["r30"] -= medians["r30"]

    del daily["month"]

    # make sure same format as reference
    daily_out["lon"]  = daily_out["lon"].astype(float)
    daily_out["lat"]  = daily_out["lat"].astype(float)
    daily_out["time"] = daily_out["time"].astype("datetime64[ns]")
    daily_out["r30"]  = daily_out["r30"].astype(np.float32)

    columns = ["lat", "lon", "time", "vx", "dx", "r30"]
    sortby = ["lat", "lon", "time"]

    daily = daily[columns].sort_values(by=sortby).reset_index(drop=True)
    daily_out = daily_out[columns].sort_values(by=sortby).reset_index(drop=True)

    pd.testing.assert_frame_equal(
        daily, daily_out, atol=1e-6
        )


    # %% visualise the distributions of the anomalies
    hist_kws = {"color": "lightgrey", "linewidth": 0.5, "edgecolor": "k",
                "bins": 50, "density": True}
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 8))
    daily.vx.plot.hist(**hist_kws, xlabel="seasonal anomaly vx (m/s)", ylabel="Density", ax=axs[0])
    daily.dx.plot.hist(**hist_kws, xlabel="seasonal anomaly dx (°)", ylabel="Density", ax=axs[1])
    daily.r30.plot.hist(**hist_kws, xlabel="seasonal anomaly r30 (m)", ylabel="Density", ax=axs[2])
    plt.tight_layout()
    
    # %% test 2: validate storm matching
    storms["time"] = pd.to_datetime(storms["date"])
    storms = storms[storms["wind"] > 0].copy()
    storms["event"] = pd.factorize(storms["cyclone_id"])[0]

    metadata = daily.merge(storms, left_on="time", right_on="time", how="inner")
    metadata = metadata[["time", "lat", "lon", "vx", "dx", "r30", "wind", "event"]]
    metadata = metadata.groupby(["lat", "lon", "event"]).max().reset_index()

    # %% how does storm max compare with recorded max winds
    metadata.plot.scatter(x="vx", y="wind", c="dx", cmap="twilight",
                          xlabel="seasonal anomaly vx (m/s)", ylabel="storm wind (m/s)",
                          title="Storm wind v ERA5 seasonal anomaly", s=1);
    if False:
        # slow (use later)
        sns.pairplot(metadata, vars=["vx", "dx", "r30", "wind"], hue="event",
                    palette="tab10", diag_kind="kde",
                    plot_kws={"alpha": 0.5, "s": 10, "edgecolor": "k"})
    # %%
    # medians  = pd.read_parquet(os.path.join(wd, medians))
    # metadata = pd.read_parquet(os.path.join(wd, metadata))
    # daily    = pd.read_parquet(os.path.join(wd, daily))
    # storms   = pd.read_parquet(os.path.join(resources_dir, storms))

    # check 0: do medians match the data_all file?
    data_all_medians = data_all.median(dim="time").to_dataframe()

    # check 1: are all storm dates present in the metadata file?

    # check 2: do winds in storms / metadata correlate with the data_all file?
    data_storms = data_all.sel(time=storms.time.values)
    data_storms = data_all.sel(time=metadata.times.values)

    # check 3: how big is the difference between storms:wind vs data_all:vx?

    # check 4: distributions of vx across all storms

    # check 5: distribution of dx across all storms

    # check 6: distributions of r30 across all storms
    
# %%