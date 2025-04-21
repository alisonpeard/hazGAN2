
rule plot_fitted_parameters:
    """
    Plot fitted parameters for each variable.
    """
    input:
        data_all=os.path.join(PROCESSING_DIR, "data_all.nc"),
        events=os.path.join(PROCESSING_DIR, "events.parquet"),
    output:
        figures=[os.path.join(FIGURE_DIR, "f01a.png")] # add b and c later
    params:
        fields=FIELDS,
        pcrit=0.05,
        cmap="PuBu_r"
    conda:
        PYENV
    script:
        os.path.join("..", "scripts", "parameter_plots.py")