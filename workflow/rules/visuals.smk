
rule plot_fitted_parameters:
    """
    Plot fitted parameters for each variable. Figure 1 in paper.
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