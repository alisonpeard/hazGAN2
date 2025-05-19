# %%
import os
import random
from environs import Env
from metrics_vis import analyze_stylegan_training
from analysis import load_samples, plot
from analysis import yflip as f

from hazGAN.plotting import figure_one, figure_two
from hazGAN.plotting import fields, samples, spatial


FIELD     = 0
MODEL     = "best"
MODEL     = str(MODEL).zfill(5) if isinstance(MODEL, int) else MODEL
THRESHOLD = 15. # None for all storms
TYPE      = "samples-1000" # ['samples', 'samples-1000]

#  begin script
if __name__ == "__main__":
    # set up environment
    env = Env()
    env.read_env()

    samples_dir = env.str("SAMPLES_DIR")
    data_dir    = env.str("DATA_DIR")

    #load training logs - - - - - - - - - - - - - - - - - - - - - - - - - - -
    try:
        metrics_path = os.path.join(samples_dir, MODEL, 'stats.jsonl')
        metrics, fig = analyze_stylegan_training(metrics_path)
    except Exception as e:
        print(e)

    # create data
    data = load_samples(samples_dir, data_dir, MODEL, threshold=THRESHOLD, sampletype=TYPE)
    u               = data['training']['uniform']
    gumbel          = data['training']['gumbel']
    x               = data['training']['data']
    mask            = data['training']['mask']

    valid_u         = data['valid']['uniform']
    valid_gumbel    = data['valid']['gumbel']
    valid_x         = data['valid']['data']
    valid_mask      = data['valid']['mask']

    samples_uniform = data['samples']['uniform']
    samples_gumbel  = data['samples']['gumbel']
    samples_x       = data['samples']['data']
    samples_mask    = data['samples']['mask']

    # %% make 64 x 64 plots of data
    if True:
        samples.plot(samples_uniform, u, field=FIELD, title="Uniform samples", cbar_label="Percentile")
        samples.plot(samples_gumbel, gumbel, field=FIELD, title="Gumbel samples")
        samples.plot(samples_x, x, field=FIELD, title="Data samples", cbar_label="Anomaly [mps]")
    
    # %% dependence between fields
    if True:
        FIELDS = [1, 2]
        fields.plot(samples_uniform, u, fields.smith1990, fields=FIELDS, figsize=.6,
                    title="Smith et al. (1990)", cbar_label="Extremal coefficient")
        fields.plot(samples_uniform, u, fields.taildependence, fields=FIELDS, figsize=.6,
                    title="Chi", cbar_label="Tail dependence", metric='chibar', thresh=.1)
        fields.plot(samples_uniform, u, fields.pearson, fields=FIELDS, figsize=.6,
                    title="Pearson correlation", cbar_label="r")

    # %% =============================================
    # %% plot correlation across space
    if False:
        # need to change these to useful scales
        # need to investigate why getting SUCH high values
        figure_two(samples_uniform, u, yflip=True, channel=0, cmap="Spectral");
        figure_two(samples_uniform, u, yflip=True, channel=1, cmap="Spectral");
        figure_two(samples_uniform, u, yflip=True, channel=2, cmap="Spectral");

        # spatial correlations
        spatial_correlations(samples_uniform, u)
    
    # %% scatter plots
    if False:
        import matplotlib.pyplot as plt
        from hazGAN.plot import plot_sample_density

        field = 0
        pixels = random.sample(range(64 * 64), 2)

        fig, axs = plt.subplots(1, 2, figsize=(6, 3))
        plot_sample_density(x[..., field], axs[0], sample_pixels=pixels, cmap="Spectral_r")
        plot_sample_density(samples_x[..., field], axs[1], sample_pixels=pixels, cmap="Spectral_r")

    # %%
