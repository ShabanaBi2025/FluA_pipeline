// modules/align_illumina.nf
process ALIGN_ILLUMINA {
    tag "$sample_id"
    label 'align'
    publishDir "${params.outdir}/aligned/illumina/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(platform), path(reads), path(fasta_ref_dir_path), val(fasta_filename), val(strain_id)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.sorted.bam"), val(strain_id), emit: bam
    tuple val(strain_id), path("${sample_id}.bwa.log"), emit: log
    tuple val(strain_id), path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    bwa mem ${fasta_ref_dir_path}/${fasta_filename} ${read1} ${read2} 2> ${sample_id}.bwa.log | samtools view -bS - | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}
