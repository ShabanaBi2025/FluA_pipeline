// modules/multiqc.nf
process MULTIQC {
    label 'multiqc_process'
    publishDir "${params.outdir}/multiqc_report", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}