"""Reference: https://github.com/alisonpeard/styleGAN-DA/blob/main/visualise.py"""
rule all_figures:
    """All figures."""
    input:
        os.path.join(FIGURE_DIR, f"{list(FIELDS.keys())[0]}.png"),
        os.path.join(FIGURE_DIR, f"{list(FIELDS.keys())[1]}.png"),
        os.path.join(FIGURE_DIR, f"{list(FIELDS.keys())[2]}.png"),
        os.path.join(FIGURE_DIR, "samples"),
        os.path.join(FIGURE_DIR, "event_intensity_barchart.png"),
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
        figa=os.path.join(FIGURE_DIR, f"{list(FIELDS.keys())[0]}.png"),
        figb=os.path.join(FIGURE_DIR, f"{list(FIELDS.keys())[1]}.png"),
        figc=os.path.join(FIGURE_DIR, f"{list(FIELDS.keys())[2]}.png")
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


rule plot_samples:
    """Figure 3: generated and observed samples."""
    input:
        train=os.path.join(TRAINING_DIR, "data.nc"),
        generated=os.path.join(GENERATED_DIR, "netcdf", "data.nc")
    output:
        outdir=directory(directory(os.path.join(FIGURE_DIR, "samples")))
    params:
        fields=FIELDS,
        shuffle=False
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
        figure=os.path.join(FIGURE_DIR, "event_intensity_barchart.png")
    params:
        month="September",
        fields=FIELDS,
        dataset=DATASET,
        do_subset=True,
        event_subset=config['event_subset']
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_barcharts.log")
    script:
        os.path.join("..", "scripts", "plot_barcharts.py")


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
        do_subset=True,
        event_subset=config['event_subset'],
        outres=16
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_correlations.log")
    script:
        os.path.join("..", "scripts", "plot_correlations.py")


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
        do_subset=True,
        event_subset=config['event_subset'],
        pois=config["points_of_interest"],
        ymin=config["latitude"]["min"],
        ymax=config["latitude"]["max"],
        xmin=config["longitude"]["min"],
        xmax=config["longitude"]["max"],
        cmap="viridis"
    conda:
        GEOENV
    log:
        file=os.path.join("logs", "plot_scatterplots.log")
    script:
        os.path.join("..", "scripts", "plot_scatterplots.py")


# rule figure_two:
#     """
#     rules:
#         - plot_fitted_parameters
#     """

# rule figure_three:
#     """
#     rules:
#         - plot_samples
#             - ws
#     """

# rule figure_four:
#     """
#     rules:
#         - plot_storm_distribution
#         - plot_spatial_correlations
#         - plot_field_correlations
#     """

# rule figure_five:
#     """
#     rules:
#         - plot_logistic_response_surface
#     """

# rule figure_six:
#     """
#     rules:
#         - plot_risk_profile
#     """

# rule figure_seven:
#     """
#     rules:
#         - plot_10_year_samples
#     """


# rule plot_storm_distribution:
#     """
#     Figure 4 (a): Saffir-Simpson storm distribution.

#     current script: /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/visualise.py
#     """

# rule plot_spatial_correlations:
#     """
#     Figure 4 (b): spatial correlations of the real and generated samples.

#     current script: /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/visualise.py
#     """

# rule plot_field_correlations:
#     """
#     Figure 4 (c): field correlations of the real and generated samples.

#     current script: /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/visualise.py
#     """

# rule plot_logistic_response_surface:
#     """
#     Figure 5: logistic response surface for each variable.

#     script: /Users/alison/Documents/DPhil/paper1.nosync/hazGAN/scripts/mangroves/train.py
#     """

# rule plot_risk_profile:
#     """
#     Figure 6: risk profile for each variable.

#     script: /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/mangrove_intersect.py
#     old: /Users/alison/Documents/DPhil/paper1.nosync/hazGAN/scripts/mangroves/riskprofile_old.py
#     """

# rule plot_10_year_samples:
#     """
#     Figure 7: 10-year samples for each variable.

#     script: /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/mangrove_intersect.py
#     """


# rule plot_brownresnick:
#     """
#     scripts:
#         - /Users/alison/Documents/DPhil/paper1.nosync/hazGAN/scripts/results/scatterplots_brownresnick.py
#         - /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/visualise.py
#     """