// modules/trimmomatic.nf
process TRIMMOMATIC {
    tag "$sample_id"
    label 'trimmomatic'
    publishDir "${params.outdir}/trimmomatic/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(platform), path(reads_list), val(strain_id)

    output:
    tuple val(sample_id), val(platform), val(strain_id), path("${sample_id}_1_paired.fq.gz"), path("${sample_id}_1_unpaired.fq.gz"), path("${sample_id}_2_paired.fq.gz"), path("${sample_id}_2_unpaired.fq.gz"), emit: reads
    tuple val(strain_id), path("*.log"), emit: log

    script:
    def read1 = reads_list[0]
    def read2 = reads_list[1]
    """
    trimmomatic PE -phred33 -summary ${sample_id}.trimmomatic_summary.log \\
        $read1 $read2 \\
        ${sample_id}_1_paired.fq.gz ${sample_id}_1_unpaired.fq.gz \\
        ${sample_id}_2_paired.fq.gz ${sample_id}_2_unpaired.fq.gz \\
        ILLUMINACLIP:${params.adapters}:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:30
    """
}