// alignment_main.nf - Main workflow for Alignment and Variant Calling

// Import common alignment processes
include { BWA_INDEX; post_alignment_processing } from './modules/align_common.nf'
// Import technology-specific alignment workflows
include { align_illumina_reads } from './modules/align_illumina.nf'
include { align_ont_reads } from './modules/align_ont.nf'


// Define channels for input cleaned reads (from the previous QC/trim pipeline)
// Illumina trimmed reads
Channel
    .fromFilePairs( params.input_trimmed_illumina, checkIfExists: true, flat: true ) { file -> file.baseName.replaceAll(/_R[12]_trimmed$/, '') }
    .ifEmpty { exit 1, "Cannot find any trimmed Illumina reads matching: ${params.input_trimmed_illumina}. Please ensure the QC/trim pipeline has been run successfully." }
    .set { trimmed_illumina_ch }

// ONT filtered reads
Channel
    .fromPath( params.input_filtered_ont, checkIfExists: true ) { file -> file.baseName.replaceAll(/_filtered\.fastq(\.gz)?$/, '') }
    .ifEmpty { exit 1, "Cannot find any filtered ONT reads matching: ${params.input_filtered_ont}. Please ensure the QC/trim pipeline has been run successfully." }
    .set { filtered_ont_ch }

// Define channels for original reference genomes (FASTA only, as GFF is not used in alignment processes)
Channel
    .fromPath( params.h1n1_reference_fasta, checkIfExists: true )
    .set { h1n1_ref_ch }

Channel
    .fromPath( params.h3n2_reference_fasta, checkIfExists: true )
    .set { h3n2_ref_ch }

workflow {
    // 1. Index H1N1 and H3N2 reference genomes (common step)
    // This will output the indexed FASTA and index files from the published directory
    all_original_references = h1n1_ref_ch.mix(h3n2_ref_ch)
    indexed_references_ch = all_original_references | BWA_INDEX

    // We need to ensure that the indexed references are available before alignment starts.
    // The `indexed_references_ch` outputs tuples like (baseName, fasta_path, index_files_glob).
    // We can use these to pass to the specific alignment modules.
    // For simplicity in calling the sub-workflows, we'll pass the original reference channels
    // and let the sub-workflows derive the indexed paths from params.output_dir.
    // This is because `BWA_INDEX` publishes, but its output channel doesn't directly
    // connect to the input of `align_illumina_reads` or `align_ont_reads` in a way
    // that carries the specific indexed file paths for each reference.
    // The `align_illumina_reads` and `align_ont_reads` workflows already derive the
    // indexed FASTA path based on sample name and params.output_dir, which is robust.

    // 2. Align Illumina reads
    aligned_illumina_bams_ch = align_illumina_reads(trimmed_illumina_ch, h1n1_ref_ch, h3n2_ref_ch)

    // 3. Align ONT reads
    aligned_ont_bams_ch = align_ont_reads(filtered_ont_ch, h1n1_ref_ch, h3n2_ref_ch)

    // 4. Combine all aligned BAMs for post-alignment processing
    // The tuple format from align_illumina_reads and align_ont_reads is (name, bam, bai, indexed_fasta)
    all_aligned_bams_ch = aligned_illumina_bams_ch.mix(aligned_ont_bams_ch)

    // 5. Run post-alignment processing (deduplication, variant calling, consensus)
    post_alignment_processing(all_aligned_bams_ch)
}
