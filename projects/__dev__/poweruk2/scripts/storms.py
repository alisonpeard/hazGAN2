"""Process storm data to create a list of storm IDs for resourcs / rfuncs (TBC)."""
# %%
import os
from glob import glob
import pandas as pd

# %%
REGION = ["Midlands", "Southwest"][0]
INDIR = "/Users/alison/Documents/DPhil/data/poweroutagesuk"

files = glob(os.path.join(INDIR, REGION, "*.txt"))
df = pd.read_csv(files[0]) # , parse_dates=["date"], index_col="date")
df.head()
# %%
# load them all into a single dataframe
dfs = []
for f in files:
    print(f)
    df = pd.read_csv(f, parse_dates=["date"], index_col="date")
    dfs.append(df)
# %%
df = pd.concat(dfs)
df = df.sort_index()
# %%
df["year"] = df.index.year
df["cyclone_id"] = df.index.year.astype(str) + "_" + df["cyclone_id"].astype(str)
df
# %%
df['cyclone_id'].nunique(), df.index.nunique()
# %%
df.to_csv(os.path.join("..", "resources", f"cyclones_{REGION.lower()}.csv"))
# %% view what the storms look like
df = pd.read_csv(os.path.join("..", "resources", f"cyclones_{REGION.lower()}.csv"), index_col="date")
df.index = pd.to_datetime(df.index)

# save these to rows of a text file
with open(os.path.join("..", "resources", f"cyclones_{REGION.lower()}.txt"), "w") as f:
    for date in df.index:
        f.write(f"{date.strftime('%Y-%m-%d')}\n")

# %%
with open(os.path.join("..", "resources", f"cyclones_{REGION.lower()}.txt"), 'r') as f:
    SELECT = [line.strip() for line in f if line.strip()]

# %%
df = df.reset_index().groupby("cyclone_number").apply(
    lambda subdf: pd.Series({
        "date": subdf.loc[subdf["wind"].idxmax(), "date"],
        "wind": subdf["wind"].max()
    })
).reset_index().set_index("date")
df = df.sort_index()
# %%
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(20, 3))

df["wind"].plot(ax=ax)
# %%
df.index

# %%
