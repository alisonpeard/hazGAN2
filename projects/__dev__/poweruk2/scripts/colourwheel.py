# %%
import matplotlib.pyplot as plt
import numpy as np

def circular_colormap(cmap='viridis'):
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    theta = np.linspace(0, 2*np.pi, 50)
    r = [0.7, 1.0]
    T, R = np.meshgrid(theta, r)
    colors = np.abs(np.sin(theta - np.pi/4))
    C, _ = np.meshgrid(colors, [0, 1])
    
    ax.pcolormesh(T, R, C, cmap=cmap)
    ax.set_rticks([])
    plt.show()

# Usage
circular_colormap('Spectral')
# %%