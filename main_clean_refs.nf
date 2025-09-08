nextflow.enable.dsl = 2

// Define references in params
params.references = [
    h1n1: file("references/h1n1/h1n1_reference.fasta", checkIfExists: true),
    h3n2: file("references/h3n2/h3n2_reference.fasta", checkIfExists: true)
]

params.outdir = params.outdir ?: "results"

workflow {
    // Convert reference map to a channel of (strain_id, fasta_path)
    Channel
        .from(params.references.collect { strain, file -> tuple(strain, file) })
        .set { original_references }

    // Clean and index reference
    CLEAN_AND_INDEX_REF(original_references)
}

process CLEAN_AND_INDEX_REF {
    tag "$strain_id"
    publishDir "${params.outdir}/cleaned_reference/${strain_id}", mode: 'copy'
    label 'align'

    input:
    tuple val(strain_id), path(original_fasta)

    output:
    tuple val(strain_id), 
          path("${strain_id}.fasta"), 
          path("${strain_id}.fasta.fai"), 
          path("${strain_id}.fasta.amb"),
          path("${strain_id}.fasta.ann"),
          path("${strain_id}.fasta.bwt"),
          path("${strain_id}.fasta.pac"),
          path("${strain_id}.fasta.sa"),
          emit: cleaned_fasta

script:
"""
set -euxo pipefail

echo "Using input FASTA: ${original_fasta}"
ls -lh "${original_fasta}"

# Clean headers and write to a new file named ${strain_id}.fasta
sed 's/ .*//' "${original_fasta}" > "${strain_id}.fasta"

# Index with samtools and BWA
samtools faidx "${strain_id}.fasta"
bwa index "${strain_id}.fasta"
"""
}
