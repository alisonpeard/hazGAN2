from pathlib import Path


rule plot_fitted_parameters:
    input:
        events=Path(PROCESSING_DIR) / "event_footprints.parquet"
    output:
        fig1=Path(FIGURE_DIR) / "parameters" / f"{list(FIELDS.keys())[0]}_upper.png",
        fig2=Path(FIGURE_DIR) / "parameters" / f"{list(FIELDS.keys())[1]}_upper.png",
        fig3=Path(FIGURE_DIR) / "parameters" / f"{list(FIELDS.keys())[2]}_upper.png"
    params:
        fields=FIELDS,
        pcrit=0.05,
        cmap="PuBu_r"
    conda:
        GEOENV
    log:
        file=Path("logs") / "plot_parameters.log"
    script:
        Path("..") / "scripts" / "plot_parameters.py"


rule plot_correlations:
    """Correlation plots of storm distribution.
    
    Swap between train=Path(TRAINING_DIR) / "data.nc" and 
    train=Path(GENERATED_DIR) / "netcdf" / "train.nc" to 
    check transforms working correctly.
    """
    input:
        train=Path(TRAINING_DIR) / "data.nc",
        # train=Path(GENERATED_DIR) / "netcdf" / "train.nc",
        generated=Path(GENERATED_DIR) / "netcdf" / "data.nc"
    output:
        dir0=directory(Path(FIGURE_DIR) / "correlations_field"),
        dir1=directory(Path(FIGURE_DIR) / "correlations_spatial")
    params:
        fields=FIELDS,
        dataset=DATASET,
        outres=8,
        subset=config['event_subset'],
        lon_min=config["longitude"]["min"],
        lon_max=config["longitude"]["max"],
        lat_min=config["latitude"]["min"],
        lat_max=config["latitude"]["max"],
        domain=config["domain"],
        plot_domain="uniform", # ["anomaly", "uniform", "standardised"]
    conda:
        GEOENV
    log:
        file=Path("logs") / "plot_correlations.log"
    script:
        Path("..") / "scripts" / "plot_correlations.py"


rule plot_samples:
    """Figure 3: generated and observed (deseasonalised) samples.
    
    Swap between train=Path(TRAINING_DIR) / "data.nc" and 
    train=Path(GENERATED_DIR) / "netcdf" / "train.nc" to 
    check transforms working correctly.
    
    >>> snakemake --profile profiles/cluster plot_samples
    """
    input:
        train=Path(TRAINING_DIR) / "data.nc",
        # train=Path(GENERATED_DIR) / "netcdf" / "train.nc",
        generated=Path(GENERATED_DIR) / "netcdf" / "data.nc"
    output:
        outdir=directory(directory(Path(FIGURE_DIR) / "samples")))
    params:
        fields=FIELDS,
        shuffle=False,
        lon_min=config["longitude"]["min"],
        lon_max=config["longitude"]["max"],
        lat_min=config["latitude"]["min"],
        lat_max=config["latitude"]["max"],
        intensity=config["event_subset"],
        domain=config["domain"],
    conda:
        GEOENV
    log:
        file=Path("logs") / "plot_samples.log"
    script:
        Path("..") / "scripts" / "plot_samples.py"


rule plot_barcharts:
    """Saffir-Simpson barcharts of storm distribution."""
    input:
        train=Path(TRAINING_DIR) / "data.nc",
        generated=Path(GENERATED_DIR) / "netcdf" / "data.nc"
    output:
        figure=Path(FIGURE_DIR) / "barcharts" / "event_intensity.png"
    params:
        event_subset=config['event_subset'],
        month="January",
        fields=FIELDS,
        dataset=DATASET
    conda:
        GEOENV
    log:
        file=Path("logs") / "plot_barcharts.log"
    script:
        Path("..") / "scripts" / "plot_barcharts.py"



rule plot_scatterplots:
    """Scatterplots of storm distribution.
    
    Need to map coords to datasets properly.
    """
    input:
        train=Path(TRAINING_DIR) / "data.nc",
        generated=Path(GENERATED_DIR) / "netcdf" / "data.nc"
    output:
        outdir=directory(Path(FIGURE_DIR) / "scatterplots")
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
        file=Path("logs") / "plot_scatterplots.log"
    script:
        Path("..") / "scripts" / "plot_scatterplots.py"


rule all_figures:
    input:
        Path(FIGURE_DIR) / "parameters" / f"{list(FIELDS.keys())[0]}_upper.png",
        Path(FIGURE_DIR) / "parameters" / f"{list(FIELDS.keys())[1]}_upper.png",
        Path(FIGURE_DIR) / "parameters" / f"{list(FIELDS.keys())[2]}_upper.png",
        Path(FIGURE_DIR) / "samples",
        Path(FIGURE_DIR) / "barcharts" / "event_intensity.png",
        Path(FIGURE_DIR) / "correlations_field",
        Path(FIGURE_DIR) / "correlations_spatial",
        Path(FIGURE_DIR) / "scatterplots"
    log:
        file=Path("logs") / "all_figures.log"