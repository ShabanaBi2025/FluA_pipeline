// modules/nextclade.nf

process EXTRACT_HA_SEGMENT {
    tag "$sample_id"
    // This is a temporary step, so we don't publish the output
    
    input:
    tuple val(sample_id), path(multi_fasta), val(strain_id)

    output:
    tuple val(sample_id), path("${sample_id}.ha.fasta"), val(strain_id)

    script:
    def ha_id = params.ha_segment_ids[strain_id]
    """
    samtools faidx ${multi_fasta}
    samtools faidx ${multi_fasta} ${ha_id} | sed '/>/! s/[RYSWKMBDHV]/N/g' > ${sample_id}.ha.fasta
    """
}

process NEXTCLADE {
    tag "$sample_id"
    label 'nextclade_process'
    publishDir "${params.outdir}/nextclade/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(ha_fasta), val(strain_id), path(dataset_dir)

    output:
    path "*.csv", emit: csv
    path "*.json", emit: json // Emit other files if needed by MultiQC or for records

    script:
    """
    nextclade run \\
        --input-dataset ${dataset_dir} \\
        --output-csv="${sample_id}.nextclade.csv" \\
        --output-json="${sample_id}.nextclade.json" \\
        ${ha_fasta}
    """
}