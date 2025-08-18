// modules/r_analysis.nf

// This process runs the R script for downstream analysis on Nextclade output.
// It requires the Nextclade CSV file and the strain ID as input.
process R_ANALYSIS {
    // Use the Conda environment defined in your new analysis.yml file.
    // This ensures that R and all required packages are available.
    conda 'envs/analysis.yml'

    // Define the inputs for the process.
    // 'strain_id' will be a string (e.g., 'h1n1', 'h3n2').
    // 'nextclade_csv' will be the path to the Nextclade output file.
    // ADDED: The path to the R script is now an input.
    input:
    tuple val(strain_id), path(nextclade_csv)
    path(r_script)

    // Define the outputs. The R script creates a summary CSV and a QC plot PNG.
    // We use dynamic naming to ensure the output files are correctly named
    // based on the strain ID.
    output:
    path("${strain_id}_clade_summary.csv"), emit: clade_summary_csv
    path("${strain_id}_qc_plot.png"), emit: qc_plot_png

    // The script block specifies the command to execute.
    // It calls your R script and passes the input file and strain ID as arguments.
    // UPDATED: The script name has been changed to 'summarise_results.R'.
    script:
    """
    Rscript "${r_script}" "${nextclade_csv}" "${strain_id}"
    """
}