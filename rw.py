# %%
import pandas as pd

path = "/Users/alison/Documents/dphil/data/ecmwf/storm-tracks/storm_track-hodges-raw-194001_202511-all-v1_0.csv"

df = pd.read_csv(path)
df["fg10"].plot.hist(bins=100)

# %%
path = "/Users/alison/Local/github/hazGAN2/projects/poweruk_winter/resources/cyclones_midlands.csv"

df = pd.read_csv(path)
df["wind"][df["wind"]>0].plot.hist(bins=100)


# %%

# %%
