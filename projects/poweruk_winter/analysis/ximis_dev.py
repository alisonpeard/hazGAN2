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


# config and paths
plt.rcParams["font.family"] = ["serif", "sans-serif", "monospace"][2]
# bd = os.path.join("/hn01-home", "spet5107") # depends on device
bd = os.path.join("/Users", "alison", "Local", "github")
wd = "hazGAN2/projects/poweruk_winter"
sys.path.append(os.path.join(bd, "hazGAN2", "workflow"))
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
    
# %% Setup
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import weibull_min, gumbel_r
from scipy.optimize import minimize_scalar

# Load your data
poi_idx = 0  # Choose a POI
field_idx = 0  # Wind speed

# Extract fitted Weibull parameters
thresh = train.isel(poi=poi_idx, field=field_idx, param=0)["params"].values.item()
scale = train.isel(poi=poi_idx, field=field_idx, param=1)["params"].values.item()
k = train.isel(poi=poi_idx, field=field_idx, param=2)["params"].values.item()

# Get exceedances (your footprint - event maxima)
x = train.isel(poi=poi_idx, field=field_idx)["anomaly"].values.flatten()
x = x[x > thresh]
x = x[~np.isnan(x)]
exceedances = x - thresh

print(f"Fitted parameters: k={k:.3f}, scale={scale:.3f}, thresh={thresh:.3f}")
print(f"Number of exceedances: {len(exceedances)}")

# %% Step 1: Transform to V^k space
# Cook & Harris (2004): "The transform of wind speed to v^k gives 
# tail-equivalence to the Exponential distribution so that extremes 
# exhibit the fastest possible convergence to the Gumbel distribution"

# %% Step 2: Calculate XIMIS plotting positions and sort data

# Sort ASCENDING (smallest first)
exceedances_sorted = np.sort(exceedances)  # Ascending: [smallest ... largest]
N = len(exceedances_sorted)
ranks = np.arange(1, N + 1)  # [1, 2, 3, ..., N]

# Now rank 1 = smallest value = smallest y
# And rank N = largest value = largest y
y_plotting = -np.log(-np.log((ranks - 0.44) / (N + 0.12)))

# Transform the sorted data
v_transformed_sorted = exceedances_sorted ** k

# %% Step 3: Fit Gumbel to transformed data
# Use the transformed sorted data for fitting

gamma_euler = 0.5772156649
sigma_gumbel = np.std(v_transformed_sorted) * np.sqrt(6) / np.pi
mu_gumbel = np.mean(v_transformed_sorted) - gamma_euler * sigma_gumbel

print(f"\nGumbel parameters (transformed space):")
print(f"  μ = {mu_gumbel:.3f}")
print(f"  σ = {sigma_gumbel:.3f}")

# %% Step 4: Visualize fit on Gumbel probability paper
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8),
                               sharex=True, sharey=True)

# Left: Original space
ax1.scatter(exceedances_sorted, y_plotting, 
            color='k', alpha=1, s=10,
            label='Data')
ax1.set_xlabel('v (original space)')
ax1.set_ylabel('Gumbel reduced variate')
ax1.set_title(f'Original space (k={k:.2f})')
ax1.grid(True, alpha=0.3)

# Right: Transformed space
ax2.scatter(v_transformed_sorted, y_plotting, 
            color='k', alpha=1, s=10,
            label='Data')
# Fitted Gumbel line
z_fit = mu_gumbel + sigma_gumbel * y_plotting
ax2.plot(z_fit, y_plotting, 'k-', linewidth=1, label='Fitted Gumbel')
ax2.set_xlabel(f'V$^{{{k:.2f}}}$ (transformed space)')
ax2.set_ylabel('Gumbel reduced variate')
ax2.set_title('Transformed space (should be more linear)')
ax2.grid(True, alpha=0.3)

for ax in (ax1, ax2):
    ax.axhline(0, color='k', linewidth=0.5, label=r"CLV $\tilde V$")
    ax.axhline(3.9, color='k', linewidth=0.5, label="1:50-yr")
    ax.axhline(9.21, color='k', linewidth=0.5, label="1:10,000-yr")
    ax.legend(loc='upper left', frameon=False)
    ax.label_outer()

plt.tight_layout()
plt.show()

# %% Step 5: Calculate return levels
# For return period RP (years), with event rate λ ≈ 1.81 events/year
# Exceedance probability: p = 1 - 1/(λ * RP)
# Gumbel reduced variate: y = -ln(-ln(p))
# Quantile in transformed space: z = μ + σ * y
# Back-transform to original space: v = z^(1/k)
# Add threshold: x = v + thresh

lambda_rate = 1.81  # Your event rate
return_periods = np.array([10, 50, 100, 500, 1000])

def cook_harris_return_level(RP, mu, sigma, k, thresh, lambda_rate):
    """
    Calculate return level using Cook & Harris methodology.
    
    From Cook et al. (2023): "For a return period RP, the exceedance 
    probability p = 1 - 1/(λ*RP) where λ is the annual event rate"
    """
    p = 1 - 1 / (lambda_rate * RP)
    y = -np.log(-np.log(p))  # Gumbel reduced variate
    z = mu + sigma * y        # Quantile in transformed space
    v = z ** (1 / k)          # Back-transform
    x = v + thresh            # Add threshold
    return x

# Calculate return levels
return_levels_new = np.array([
    cook_harris_return_level(rp, mu_gumbel, sigma_gumbel, k, thresh, lambda_rate)
    for rp in return_periods
])

# Compare with your current Weibull approach
def weibull_return_level(RP, k, scale, thresh, lambda_rate):
    """Your current approach: direct Weibull quantiles"""
    p = 1 - 1 / (lambda_rate * RP)
    v = scale * (-np.log(1 - p)) ** (1 / k)  # Weibull quantile
    return v + thresh

return_levels_old = np.array([
    weibull_return_level(rp, k, scale, thresh, lambda_rate)
    for rp in return_periods
])

# %% Step 6: Compare approaches
print("\nReturn Level Comparison:")
print(f"{'RP (years)':<12} {'Weibull':<12} {'Cook-Harris':<12} {'Difference %':<12}")
print("-" * 50)
for i, rp in enumerate(return_periods):
    diff_pct = 100 * (return_levels_new[i] - return_levels_old[i]) / return_levels_old[i]
    print(f"{rp:<12} {return_levels_old[i]:<12.2f} {return_levels_new[i]:<12.2f} {diff_pct:+.2f}%")

# %% Step 7: Implement for your invPIT function
def gumbel_to_original(u, mu, sigma, k, thresh):
    """
    Transform uniform [0,1] to original space via Gumbel + back-transform.
    This replaces your Weibull quantile function.
    """
    y = -np.log(-np.log(u))  # Uniform to Gumbel reduced variate
    z = mu + sigma * y        # Gumbel quantile in transformed space
    v = z ** (1 / k)          # Back-transform
    x = v + thresh            # Add threshold
    return x

# Test round-trip
u_test = train["uniform"].isel(poi=poi_idx, field=field_idx).values
x_reconstructed = gumbel_to_original(u_test, mu_gumbel, sigma_gumbel, k, thresh)

# Plot comparison
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(x[x > thresh], bins=30, alpha=0.5, density=True, label='Original exceedances')
ax.hist(x_reconstructed[x_reconstructed > thresh], bins=30, alpha=0.5, 
        density=True, label='Reconstructed (Cook-Harris)')
ax.axvline(thresh, color='k', linestyle='--', alpha=0.3)
ax.set_xlabel('Wind speed')
ax.set_ylabel('Density')
ax.legend()
ax.set_title(f'Round-trip test for {train["poi"].values[poi_idx].title()}')
plt.show()
# %%
