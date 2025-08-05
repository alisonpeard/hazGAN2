"""Reference: https://github.com/alisonpeard/styleGAN-DA/blob/main/visualise.py"""
rule all_figures:
    """All figures."""
    input:
        os.path.join(FIGURE_DIR, "parameters", f"{list(FIELDS.keys())[0]}_upper.png"),
        os.path.join(FIGURE_DIR, "parameters", f"{list(FIELDS.keys())[1]}_upper.png"),
        os.path.join(FIGURE_DIR, "parameters", f"{list(FIELDS.keys())[2]}_upper.png"),
        os.path.join(FIGURE_DIR, "samples"),
        os.path.join(FIGURE_DIR, "barcharts", "event_intensity.png"),
        os.path.join(FIGURE_DIR, "correlations_field"),
        os.path.join(FIGURE_DIR, "correlations_spatial"),
        os.path.join(FIGURE_DIR, "scatterplots")
    log:
        file=os.path.join("logs", "all_figures.log")


rule plot_fitted_parameters:
    """
    Figure 2: fitted parameters for each variable. Figure 1 in paper.

    >>> snakemake --profile profiles/cluster plot_fitted_parameters
    """
    input:
        events=os.path.join(PROCESSING_DIR, "events.parquet")
    output:
        fig1=os.path.join(FIGURE_DIR, "parameters", f"{list(FIELDS.keys())[0]}_upper.png"),
        fig2=os.path.join(FIGURE_DIR, "parameters", f"{list(FIELDS.keys())[1]}_upper.png"),
        fig3=os.path.join(FIGURE_DIR, "parameters", f"{list(FIELDS.keys())[2]}_upper.png")
    params:
        fields=FIELDS,
        pcrit=0.05,
        cmap="PuBu_r"
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_parameters.log")
    script:
        os.path.join("..", "scripts", "plot_parameters.py")


rule plot_correlations:
    """Correlation plots of storm distribution."""
    input:
        train=os.path.join(TRAINING_DIR, "data.nc"),
        generated=os.path.join(GENERATED_DIR, "netcdf", "data.nc")
    output:
        dir0=directory(os.path.join(FIGURE_DIR, "correlations_field")),
        dir1=directory(os.path.join(FIGURE_DIR, "correlations_spatial"))
    params:
        fields=FIELDS,
        dataset=DATASET,
        outres=32,
        event_subset=config['event_subset'],
        lon_min=config["longitude"]["min"],
        lon_max=config["longitude"]["max"],
        lat_min=config["latitude"]["min"],
        lat_max=config["latitude"]["max"],
        domain="standardised", # ["anomaly", "uniform", "standardised"], "uniform needed for Smith (1990)"
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_correlations.log")
    script:
        os.path.join("..", "scripts", "plot_correlations.py")


rule plot_samples:
    """Figure 3: generated and observed samples.
    
    >>> snakemake --profile profiles/cluster plot_samples
    """
    input:
        train=os.path.join(TRAINING_DIR, "data.nc"),
        generated=os.path.join(GENERATED_DIR, "netcdf", "data.nc")
    output:
        outdir=directory(directory(os.path.join(FIGURE_DIR, "samples")))
    params:
        fields=FIELDS,
        shuffle=False,
        lon_min=config["longitude"]["min"],
        lon_max=config["longitude"]["max"],
        lat_min=config["latitude"]["min"],
        lat_max=config["latitude"]["max"],
        domain=config["domain"],
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_samples.log")
    script:
        os.path.join("..", "scripts", "plot_samples.py")


rule plot_barcharts:
    """Saffir–Simpson barcharts of storm distribution."""
    input:
        train=os.path.join(TRAINING_DIR, "data.nc"),
        generated=os.path.join(GENERATED_DIR, "netcdf", "data.nc")
    output:
        figure=os.path.join(FIGURE_DIR, "barcharts", "event_intensity.png")
    params:
        month="September",
        fields=FIELDS,
        dataset=DATASET,
        event_subset=config['event_subset']
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_barcharts.log")
    script:
        os.path.join("..", "scripts", "plot_barcharts.py")



rule plot_scatterplots:
    """Scatterplots of storm distribution.
    
    Need to map coords to datasets properly...
    """
    input:
        train=os.path.join(TRAINING_DIR, "data.nc"),
        generated=os.path.join(GENERATED_DIR, "netcdf", "data.nc")
    output:
        outdir=directory(os.path.join(FIGURE_DIR, "scatterplots"))
    params:
        event_subset=config['event_subset'],
        pois=config["points_of_interest"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        cmap="viridis",
        channel_labels={field: FIELDS[field]["title"] for field in FIELDS.keys()}
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_scatterplots.log")
    script:
        os.path.join("..", "scripts", "plot_scatterplots.py")