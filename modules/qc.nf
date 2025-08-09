// modules/qc.nf
process FASTQC {
    tag "$sample_id"
    label 'fastqc'
    publishDir "${params.outdir}/fastqc/${strain_id}/${platform}", mode: 'copy'

    input:
    tuple val(sample_id), val(platform), path(reads), val(strain_id)

    output:
    path("*.{zip,html}"), emit: reports
    
    script:
    def fastqc_inputs = (reads instanceof List) ? reads.join(' ') : reads
    """
    fastqc $fastqc_inputs -o ./
    """
}