# %%
import os
import numpy as np
import pandas as pd
import xarray as xr
from calendar import month_name
import matplotlib.pyplot as plt


wd = os.path.join("..", "..", "results", "processing")
resources_dir = os.path.join("..", "..", "resources")

metadata = "event_metadata.parquet"
fitted   = "events.parquet"


if __name__ == "__main__":
    # load the data
    metadata = pd.read_parquet(os.path.join(wd, metadata))
    fitted = pd.read_parquet(os.path.join(wd, fitted))

0# %%