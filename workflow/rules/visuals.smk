rule plot_fitted_parameters:
    """
    Figure 2: fitted parameters for each variable. Figure 1 in paper.
    """
    input:
        data_all=os.path.join(PROCESSING_DIR, "data_all.nc"),
        events=os.path.join(PROCESSING_DIR, "events.parquet"),
    output:
        figa=os.path.join(FIGURE_DIR, "f01a.png"),
        figb=os.path.join(FIGURE_DIR, "f01b.png"),
        figc=os.path.join(FIGURE_DIR, "f01c.png")
    params:
        fields=FIELDS,
        pcrit=0.05,
        cmap="PuBu_r"
    conda:
        PYENV
    log:
        file=os.path.join("logs", "plot_parameters.log")
    script:
        os.path.join("..", "scripts", "parameter_plots.py")


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
# rule plot_samples:
#     """
#     Figure 3: generated and observed samples.

#     wildcards:
#         - field: {field}

#     current script: /Users/alison/Documents/DPhil/github.nosync/styleGAN-DA/visualise.py
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