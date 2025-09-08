// modules/snpeff.nf
process RUN_SNP_EFF {
    label 'snpeff'
    conda "${baseDir}/envs/environment.yml"

    input:
    tuple val(sample_id), path(vcf_file), val(genome), val(subtype)

    output:
    path "${sample_id}.annotated.vcf"
    path "snpEff_genes.txt"
    path "snpEff_summary.html"

    publishDir "results/snpeff/output/${subtype}/${sample_id}", mode: 'copy'

    script:
    """
    snpEff -c ${baseDir}/results/snpeff/data/snpEff.config $genome $vcf_file > ${sample_id}.annotated.vcf
    """
}
