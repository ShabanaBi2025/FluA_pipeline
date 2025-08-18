// modules/bcftools.nf
process BCFTOOLS_MPILEUP {
    tag "$sample_id"
    label 'align'
    publishDir "${params.outdir}/variants/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(reference_fasta), val(strain_id)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi"), path(reference_fasta), val(strain_id), emit: vcf
    tuple val(strain_id), path("${sample_id}.vcf.stats"), emit: stats

    script:
    """
    # Indel-aware mpileup with allele depth annotations
    bcftools mpileup -Ou -f ${reference_fasta} -a AD,DP,SP ${sorted_bam} | \
        bcftools call -mv --ploidy 1 -Oz -o ${sample_id}.vcf.gz

    # Index and stats
    bcftools index ${sample_id}.vcf.gz
    bcftools stats ${sample_id}.vcf.gz > ${sample_id}.vcf.stats
    """
}

process BCFTOOLS_CONSENSUS {
    tag "$sample_id"
    label 'align'
    publishDir "${params.outdir}/consensus/${strain_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_gz), path(vcf_index), path(reference_fasta), val(strain_id)

    output:
    tuple val(sample_id), path("${sample_id}.consensus.fasta"), val(strain_id), emit: fasta

    script:
    """
    if [ \$(bcftools view -H ${vcf_gz} | wc -l) -eq 0 ]; then
        cp ${reference_fasta} ${sample_id}.consensus.fasta
    else
        bcftools consensus -f ${reference_fasta} ${vcf_gz} -o ${sample_id}.consensus.fasta
    fi
    """
}