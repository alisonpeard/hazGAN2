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


if __name__ == "__main__":
    # load the data
    data_all = xr.open_dataset(os.path.join(wd, data_all))
    storms = pd.read_csv(os.path.join(resources_dir, storms))
    metadata = pd.read_parquet(os.path.join(wd, metadata))
    medians_out = pd.read_parquet(os.path.join(wd, medians))
    daily_out = pd.read_parquet(os.path.join(wd, daily))

    # %% check 1: validate deseasonalisation
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

    # %% 
    # PLAN:
    # (a) fit distributions to both min/max tails 
    # (b) transform uniform to double-gumbel
    hist_kws = {"linewidth": 0.5, "edgecolor": "k",
                "bins": 50, "density": True}

    u = np.array(sorted(np.random.uniform(size=10_000)))

    def laplace_quantile(u, mu=0, b=1):
        return mu + b * np.sign(u - 0.5) * np.log(1 - 2 * np.abs(u - 0.5))

    def gumbel_quantile(u, mu=1, beta=1):
        return mu - beta * np.log(-np.log(u))
    
    shift = 0
    y0 = gumbel_quantile(u, mu=shift)
    y1 = -gumbel_quantile(1 - u, mu=shift)
    y3 = 0.5 * (y0 + y1)
    y = laplace_quantile(u)

    fig, axs = plt.subplots(2, 3, figsize=(15, 8))

    ax = axs[0, 0]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y0, **hist_kws, label="Gumbel", color="red", alpha=0.5);
    ax.hist(y1, **hist_kws, label="Gumbel (negative)", color="yellow", alpha=0.5);
    ax.set_xlim(-10, -5)
    ax.set_ylim(0, 0.01)

    ax = axs[0, 1]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y0, **hist_kws, label="Gumbel", color="red", alpha=0.5);
    ax.hist(y1, **hist_kws, label="Gumbel (negative)", color="yellow", alpha=0.5);
    # ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="blue", alpha=0.5);
    ax = axs[0, 2]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y0, **hist_kws, label="Gumbel", color="red", alpha=0.5);
    ax.hist(y1, **hist_kws, label="Gumbel (negative)", color="yellow", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="blue", alpha=0.5);
    ax.set_xlim(5, 10)
    ax.set_ylim(0, 0.01)
    ax.legend()

    ax = axs[1, 0]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="blue", alpha=0.5);
    ax.set_xlim(-10, -5)
    ax.set_ylim(0, 0.01)

    ax = axs[1, 1]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="blue", alpha=0.5);

    ax = axs[1, 2]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="blue", alpha=0.5);
    ax.set_xlim(5, 10)
    ax.set_ylim(0, 0.01)
    ax.legend()
    # %%
    # check 2: validate storm matching

    # check 3: how big is the difference between storms:wind vs data_all:u10?

    # check 4: distributions of u10 across all storms

    # check 5: distribution of v10 across all storms

    # check 6: distributions of r30 across all storms
    
# %%