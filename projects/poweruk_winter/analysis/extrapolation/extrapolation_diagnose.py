# %%
import os
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

os.chdir('/hn01-home/spet5107/hazGAN2')


def gumbel(uniform):
    """uniform -> Gumbel(0, 1)"""
    maxval = np.max(uniform)
    if maxval > 1:
        raise ValueError(f"Some uniform > 1 ({maxval})")
    return -np.log(-np.log(uniform))


def laplace(uniform, mu=0, b=1):
    """uniform -> Laplace(mu, b) (quantile function)"""
    maxval = np.max(uniform)
    if maxval > 1:
        raise ValueError(f"Some uniform > 1 ({maxval})")
    
    return np.where(
        uniform <= 0.5, 
        mu + b * np.log(2 * uniform),
        mu - b * np.log(2 - 2 * uniform)
        )

def inv_gumbel(x):
    """Gumbel(0, 1) -> uniform"""
    return np.exp(-np.exp(-x))

def inv_laplace(x, mu=0, b=1):
    """Laplace(mu, b) -> uniform (CDF function)."""
    return np.where(
        x <= mu,
        0.5 * np.exp((x - mu) / b),
        1 - 0.5 * np.exp(-(x - mu) / b)
    )

k = 1
project = ["bayofbengal_era5", "poweruk_winter"][k]
inv_transform = [inv_gumbel, inv_laplace][k]

trn_path = os.path.join("projects", project, "results", "training", "data.nc")
stat_dir = os.path.join("projects", project, "results", "training", "image_stats.npz")
trn_dir = os.path.join("projects", project, "results", "training", "rgb")
gen_dir = os.path.join("projects", project, "results", "generated", "images")

os.listdir(gen_dir)
# %%
import xarray as xr

rp_max = 10_000
p_min = 1 / rp_max
p_max = 1 - 1 / rp_max

train = xr.open_dataset(trn_path)
u = train.uniform.values

assert (u < p_max).all(), "Found value above {}-year return level".format(rp_max)
assert (p_min < u).all(), "Found value below {}-year return level".format(rp_max)

y = laplace(u)
y_lower = laplace(p_min)
y_upper = laplace(p_max)
z = (y - y_lower) / (y_upper - y_lower)

# %%



# %%
print("{} contains {} samples".format(trn_dir, len(os.listdir(trn_dir))))
print("{} contains {} samples".format(gen_dir, len(os.listdir(gen_dir))))
# %%

def load_npys(dir):
    samples = []
    for fname in tqdm(os.listdir(dir)):
        if fname.endswith(".npy") is False:
            continue
        sample = np.load(os.path.join(dir, fname))
        samples.append(sample / 255.0)
    samples = np.stack(samples, axis=0)
    return samples

def load_pngs(dir):
    from PIL import Image
    samples = []
    for fname in tqdm(os.listdir(dir)):
        if fname.endswith(".png") is False:
            continue
        sample = np.array(Image.open(os.path.join(dir, fname)))
        samples.append(sample / 255.0)
    samples = np.stack(samples, axis=0)
    return samples

if os.listdir(gen_dir)[0].endswith(".npy"):
    print("Loading npys")
    load_fn = load_npys
elif os.listdir(gen_dir)[0].endswith(".png"):
    print("Loading pngs")
    load_fn = load_pngs

train = load_fn(trn_dir)
gen = load_fn(gen_dir)

# %%
# apply image statistics to rescale
image_stats = np.load(stat_dir)
image_minima = image_stats['min']
image_maxima = image_stats['max']
image_n      = image_stats['n']
image_range  = image_maxima - image_minima

# rescale images to marginals scale
train_y = (train * (image_n + 1) - 1) / (image_n - 1)
train_y = train_y * image_range + image_minima
train_u = inv_transform(train_y)

gen_y = (gen * (image_n + 1) - 1) / (image_n - 1)
gen_y = gen_y * image_range + image_minima
gen_u = inv_transform(gen_y)
# %%
train_maxima = np.max(train, axis=(0))
gen_maxima = np.max(gen, axis=(0))
# %%
fig, ax = plt.subplots(figsize=(10, 5))
ax.hist(train_y[..., 0].ravel(), density=True, alpha=0.5, label="train_y")
ax.hist(gen_y[..., 0].ravel(), density=True, alpha=0.5, label="gen_y")
ax.legend(loc="upper right")
# %%
print("Scaled: {}, {}".format(train[..., 0].max(), gen[..., 0].max()))
print("Laplace: {}, {}".format(train_y[..., 0].max(), gen_y[..., 0].max()))
print("Uniform: {}, {}".format(train_u[..., 0].max(), gen_u[..., 0].max()))

# %%
print("Scaled: {}, {}".format(train[..., 0].min(), gen[..., 0].min()))
print("Laplce: {}, {}".format(train_y[..., 0].min(), gen_y[..., 0].min()))
print("Uniform: {}, {}".format(train_u[..., 0].min(), gen_u[..., 0].min()))
# %%

# %%


images_x = stats.invPIT(images_u, train_x, theta=theta, distns=distns, two_tailed=two_tailed)