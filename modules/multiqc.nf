// modules/multiqc.nf
process MULTIQC {
    tag "$strain_id" // Tag the process with the strain name
    label 'multiqc_process'
    publishDir "${params.outdir}/multiqc_report/${strain_id}", mode: 'copy'

    input:
    tuple val(strain_id), path(files) // Input is now a tuple: (strain_id, list_of_files)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}