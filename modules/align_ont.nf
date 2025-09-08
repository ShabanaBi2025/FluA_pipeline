// modules/align_ont.nf
process ALIGN_ONT {
    tag "$sample_id"
    label 'align'
    publishDir "${params.outdir}/aligned/ont/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(platform), path(filtered_read), path(fasta_ref_dir_path), val(fasta_filename), val(strain_id)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.sorted.bam"), val(strain_id), emit: bam
    tuple val(sample_id), val(platform), path("${sample_id}.sorted.bam.bai"), val(strain_id), emit: bai
    tuple val(strain_id), path("${sample_id}.minimap2.log"), emit: log
    tuple val(strain_id), path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont --MD \\
        ${fasta_ref_dir_path}/${fasta_filename} ${filtered_read} \\
        2> ${sample_id}.minimap2.log | \\
        samtools view -bS - | \\
        samtools sort -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}