// modules/nextclade.nf
process EXTRACT_HA_SEGMENT {
    tag "$sample_id"

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

process NEXTCLADE_STRAIN {
    tag "$strain_id"
    label 'nextclade'
    publishDir "${params.outdir}/nextclade/${strain_id}", mode: 'copy'

    input:
    tuple val(strain_id), path(ha_fastas), path(dataset_dir)

    output:
    tuple val(strain_id), path("${strain_id}_nextclade_results"), emit: results

    script:
    """
    # Create a directory for the combined results
    mkdir -p ${strain_id}_nextclade_results

    # Concatenate all individual HA fasta files into one file
    cat ${ha_fastas.join(' ')} > combined_samples.fasta

    # Run Nextclade on the combined file
    nextclade run \\
      --input-dataset ${dataset_dir} \\
      --output-all=${strain_id}_nextclade_results/ \\
      combined_samples.fasta
    """
}
