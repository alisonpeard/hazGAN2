"""
Dev for two-tailed extremes.
"""
# %%
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    # Testing training distribution
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

    ax = axs[0, 2]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y0, **hist_kws, label="Gumbel", color="red", alpha=0.5);
    ax.hist(y1, **hist_kws, label="Gumbel (negative)", color="yellow", alpha=0.5);
    ax.set_xlim(5, 10)
    ax.set_ylim(0, 0.01)
    ax.legend()

    ax = axs[1, 0]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="brown", alpha=0.5);
    ax.set_xlim(-10, -5)
    ax.set_ylim(0, 0.01)

    ax = axs[1, 1]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="brown", alpha=0.5);

    ax = axs[1, 2]
    ax.hist(y, **hist_kws, label="Laplace", color="orange", alpha=0.5);
    ax.hist(y3, **hist_kws, label="Gumbel (difference)", color="brown", alpha=0.5);
    ax.set_xlim(5, 10)
    ax.set_ylim(0, 0.01)
    ax.legend()

# %%