# %%
import os
import sys
import yaml
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def format_poi(p):
    """Handle different formats of points of interest"""
    if isinstance(p, list):
        p = {"lon": p[0], "lat": p[1]}
    return p


def format_pois(p: list):
    """Format for extraction from xarray"""
    keys = list(p.keys())
    p = [format_poi(p_i) for p_i in p.values()]
    lats = [p_i["lat"] for p_i in p]
    lons = [p_i["lon"] for p_i in p]
    print(f"{lats=}, {lons=}")
    return keys, lats, lons


if __name__ == "__main__":
    # config and paths
    plt.rcParams["font.family"] = ["serif", "sans-serif", "monospace"][2]
    bd = os.path.join("/hn01-home", "spet5107") # depends on device
    wd = "hazGAN2/projects/poweruk_winter"
    sys.path.append(os.path.join(bd, "hazGAN2", "workflow"))
    from src.python.statistics import invPIT
    os.chdir(bd)

    with open(os.path.join(wd, "config.yaml"), "r") as stream:
        config = yaml.safe_load(stream)

    fields = list(config["fields"].keys())

    train = xr.open_dataset(os.path.join(wd, "results", "testing", "pois.nc"))
    pois = train["poi"].values
    num_pois = len(pois)
    print("Loaded pois from file: {}".format(os.path.join(wd, "results", "testing", "pois.nc")))
    print("Pois: {}".format(pois))
    print("Train for poi {}:\n\n{}".format(pois[0].title(), train.sel(poi=pois[0])))
    # %% view distributions overall
    k = 0

    fig, axs = plt.subplots(1, num_pois, figsize=(3*num_pois, 2), constrained_layout=True)

    for i, poi in enumerate(pois):
        ax = axs[i] if num_pois > 1 else axs
        data = train.isel(
            field=k, poi=i
            )["anomaly"].values.flatten()
        data = data[~np.isnan(data)]
        ax.hist(data, bins=30, density=True, alpha=0.5, label=fields[k])
        ax.set_title("{} (n={})".format(poi.title(), len(data)))
        ax.legend(loc="upper left", frameon=False)
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")

    full_hist = (fig, axs)

    # %% view exceedances

    fig, axs = plt.subplots(1, num_pois, figsize=(3*num_pois, 2), constrained_layout=True)

    for i, poi in enumerate(pois):
        ax = axs[i] if num_pois > 1 else axs

        loc = data = train.isel(
            field=k, poi=i, param=0
            )["params"].values.item()
        
        data = train.isel(
            field=k, poi=i
            )["anomaly"].values.flatten()
        
        data = data[data > loc]
        data = data[~np.isnan(data)]
        ax.hist(data, bins=30, density=True, alpha=0.5, label=fields[k])
        ax.set_title("{} (n={})".format(poi.title(), len(data)))
        ax.legend(loc="upper left", frameon=False)
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")

    tail_hist = (fig, axs)

    # %% do invPIT from uniform and Laplce

    u = train["uniform"].values
    x = train["anomaly"].values
    theta = train["params"].values

    domain = config["domain"]
    distns = [v["distn"] for v in config["fields"].values()]
    two_tailed = [v["two_tailed"] for v in config["fields"].values()]

    print("array shapes: u={}, x={}, θ={}".format(u.shape, x.shape, theta.shape))
    print(f"{domain=}, {distns=}, {two_tailed=}")

    # %% invPIT

    x_inv = invPIT(
        u, x, theta=theta, domain=None, distns=distns, two_tailed=two_tailed
        )
    print("invPIT shape: {}".format(x_inv.shape))
    # %%
    axs = full_hist[1]

    for i, poi in enumerate(pois):
        ax = axs[i] if num_pois > 1 else axs
        data = x_inv[:, i, k]
        data = data[~np.isnan(data)]
        ax.hist(data, bins=30, density=True, alpha=0.5, label=r"$T^{{-1}}(T(x))$")
        ax.set_title("{} (n={})".format(poi.title(), len(data)))
        ax.legend(loc="upper left", frameon=False)
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")

    full_hist[0]

    # %% 
    axs = tail_hist[1]

    for i, poi in enumerate(pois):
        ax = axs[i] if num_pois > 1 else axs

        loc = train.isel(
            field=k, poi=i, param=0
            )["params"].values.item()
        
        data = x_inv[:, i, k]
        
        data = data[data > loc]
        data = data[~np.isnan(data)]
        ax.hist(data, bins=30, density=True, alpha=0.5, label=r"$T^{{-1}}(T(x))$")
        ax.set_title("{} (n={})".format(poi.title(), len(data)))
        ax.legend(loc="upper left", frameon=False)
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")

    tail_hist[0]

    # %%
    rp_max = 100_000
    p_max = 1 - 1/rp_max

    u_max = np.full(u.shape, p_max)
    x_max = invPIT(
        u_max, x, theta=theta, domain=None, distns=distns, two_tailed=two_tailed
        )
    print("invPIT max shape: {}".format(x_max.shape))
    # %%    
    axs = tail_hist[1]

    for i, poi in enumerate(pois):
        ax = axs[i] if num_pois > 1 else axs

        loc = train.isel(
            field=k, poi=i, param=0
            )["params"].values.item()
        
        data = x_max[:, i, k]
        
        data = data[data > loc]
        data = data[~np.isnan(data)]
        ax.axvline(data[0], color="r", linestyle="--", label=rf"$T^{{-1}}$({rp_max:,.0f}-yr)")
        ax.set_title("{} (n={})".format(poi.title(), len(data)))
        ax.legend(loc="upper left", frameon=False)
        ax.set_xlabel("Value")
        ax.set_ylabel("Density")
    
    tail_hist[0]

# %% view transforms
