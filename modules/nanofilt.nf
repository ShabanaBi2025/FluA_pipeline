// modules/nanofilt.nf
process NANOFILT {
    tag "$sample_id"
    label 'nanofilt'
    publishDir "${params.outdir}/nanofilt/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(platform), path(read), val(strain_id)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}_filtered.fastq.gz"), val(strain_id), emit: reads

    script:
    """
    zcat ${read} | NanoFilt -q 10 -l 500 | gzip > ${sample_id}_filtered.fastq.gz
    """
}
