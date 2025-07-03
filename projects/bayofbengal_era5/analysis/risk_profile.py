"""
Plot the Lamb (2010) risk curve for mangrove damages.
"""
# %%
import os
import xarray as xr
import matplotlib.pyplot as plt

plt.style.use('default')
plt.rcParams['legend.frameon'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'

blues = ['#fff2ccff', '#a2c4c9ff', '#0097a7ff', '#0b5394ff']


def riskprofileplot(
        tree:xr.DataTree, label:str, ax:plt.Axes, truncate:bool=False,
        minrp:float=1, maxrp:float=500, **kwargs
        ):
    ds = tree[label].to_dataset()
    ds = ds.where(ds['return_period'] >= minrp, drop=True)
    ds = ds.where(ds['return_period'] <= maxrp, drop=True)
    ax.plot(ds['return_period'], ds['expected_damage'], label=label, **kwargs)


if __name__ == "__main__":
    
    wd = os.path.join("..", "results", "mangroves")
    # tree = xr.open_dataset(os.path.join(wd, "damage_scenarios.nc"))#
    # tree = xr.DataTree(tree)

    tree = xr.open_datatree(os.path.join(wd, "damage_scenarios.nc"))

    # plot parameters
    scatter_kwargs = {
        'linestyle': 'none', 'marker': 'o', 'mfc': 'k', 'mec': 'k',
        'mew': 0.25, 'alpha': 0.8, 'ms': 4
        }
    line_kwargs = {
        'linewidth': 2, 'alpha': 0.8, 'marker': 'o', 'mfc': 'none', 'ms': 4
        }

    xmin = min([ds['return_period'].min() for ds in tree.values()]).data.item()
    xmax = max([ds['return_period'].max() for ds in tree.values()]).data.item()
    ymin = min([ds['expected_damage'].min() for ds in tree.values()]).data.item()
    ymax = max([ds['expected_damage'].max() for ds in tree.values()]).data.item()


    # plot the risk curve
    fig, ax = plt.subplots(figsize=(7.0, 4.33))

    riskprofileplot(tree, 'ERA5', ax, color='k', **scatter_kwargs)
    riskprofileplot(tree, 'HazGAN', ax, color="#4682B4", **line_kwargs)
    riskprofileplot(tree, 'Independent', ax, color=blues[3], zorder=0, **line_kwargs)
    riskprofileplot(tree, 'Dependent', ax, color=blues[2], **line_kwargs)


    # legend options
    # ax.legend(bbox_to_anchor=(1.05, .8), loc='upper left',
    #           title="Dataset")  # Places legend outside
    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),  # Center horizontally, below the plot
        ncol=4,  # Spread entries horizontally
        frameon=False,
        handletextpad=0.5,  # Reduce space between handle and text (default is 0.8)
        columnspacing=2.0,  # Reduce space between columns (default is 2.0)
        labelspacing=0.4
        # title="Dataset"
    )
    plt.setp(ax.get_legend().get_title(), fontsize='12', fontweight='bold')

    # def annotate(line, label, above=True):
    #     x = line.get_xdata()
    #     y = line.get_ydata()[-1]
    #     yloc = 5 if above else -20
    #     plt.annotate(label, xy=(x[-1], y), xytext=(-20, yloc), 
    #                  textcoords='offset points', va='center', ha='right',
    #                  color=line.get_color(), fontsize=12, fontweight='bold')

    # annotate(plt.gca().lines[0], 'ERA5')
    # annotate(plt.gca().lines[1], 'HazGAN')
    # annotate(plt.gca().lines[2], 'Independent')
    # annotate(plt.gca().lines[3], 'Dependent', above=False)

    # configure x-axis
    ax.set_xscale('log')
    yticks = [ymin, 2000, 4000, 6000, ymax]
    yticklabels = ["", "2000", "4000", "6000", f"({ymax:.0f})"]

    ax.set_yticks(yticks, labels=yticklabels)
    ax.spines['left'].set_bounds(ymin, ymax)
    ax.text(0, ymin - 400, f"({ymin:.0f})", transform=ax.get_yaxis_transform(), 
            ha='right', va='center', fontsize=12)

    ax.set_xticks([1, 10, 100, 500], labels=['1', '10', '100', '500'])
    ax.spines['bottom'].set_bounds(1, 500)

    # ticks config
    ax.tick_params(direction='in')
    ax.xaxis.set_tick_params(which='minor', direction='in')
    ax.yaxis.set_tick_params(which='minor', direction='in')
    ax.tick_params(axis='x', length=8)
    ax.tick_params(axis='y', length=8)
    ax.tick_params(axis='x', which='minor', length=4)
    ax.tick_params(axis='y', which='minor', length=4)

    # turn off minor x-ticks
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.xaxis.set_minor_locator(plt.NullLocator())

    # mark the 1-year return period
    # ax.fill_betweenx([ymin, ymax], 0, 1, color="#F1F3F5", alpha=0.8, zorder=0) #'#F4F1EA'
    # ax.axvline(x=1, ymax=0.95, color='#333333', linestyle='dashed', linewidth=1, zorder=1)

    ax.set_xlabel("Return period (years)", fontsize=13, fontweight='bold')
    # ax.set_ylabel("Expected damage area (km$^2$)", fontsize=13, fontweight='bold')
    ax.set_ylabel("Expected\ndamage\narea\n(km$^2$)", 
                fontsize=12, 
                fontweight='bold',
                rotation=0,
                labelpad=15,  # Reduced padding
                va='center',
                ha='right')

    plt.subplots_adjust(left=0.2) 
    plt.tight_layout()

    plt.tight_layout()  
    os.makedirs('../results/figures/mangroves', exist_ok=True)
    plt.savefig('../results/figures/mangroves/risk_profile.pdf', dpi=300, bbox_inches='tight')

# %%