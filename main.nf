// main.nf -  Workflow

nextflow.enable.dsl = 2

// === INCLUDE ALL PROCESS MODULES ===
include { FASTQC }             from './modules/qc.nf'
include { TRIMMOMATIC }        from './modules/trimmomatic.nf'
include { NANOFILT }           from './modules/nanofilt.nf'
include { ALIGN_ILLUMINA }     from './modules/align_illumina.nf'
include { ALIGN_ONT }          from './modules/align_ont.nf'
include { BCFTOOLS_MPILEUP }   from './modules/bcftools.nf'
include { BCFTOOLS_CONSENSUS } from './modules/bcftools.nf'
include { EXTRACT_HA_SEGMENT } from './modules/nextclade.nf'
include { NEXTCLADE }          from './modules/nextclade.nf'
include { MULTIQC }            from './modules/multiqc.nf'

// === FULL WORKFLOW ===
workflow {

    // --- 1. QC and Trimming ---
    Channel
        .fromPath(params.raw_reads)
        .map { file ->
            def strain_id = file.parent.name
            def sample_id = file.baseName.replaceAll(/_R?[12]|\.fastq(.gz)?/, '')
            def platform = (file.name.contains('_R1') || file.name.contains('_R2')) ? 'illumina' : 'ont'
            return tuple(strain_id, sample_id, platform, file)
        }
        .groupTuple(by: [0, 1])
        .map { strain_id, sample_id, platforms, files ->
            def platform = platforms[0]
            def reads = (platform == 'illumina') ? files.sort() : files[0]
            return tuple(sample_id, platform, reads, strain_id)
        }
        .set { reads_ch }

    FASTQC(reads_ch)

    reads_ch
        .branch { sample_id, platform, reads, strain_id ->
            illumina: platform == 'illumina'
            ont:      platform == 'ont'
        }
        .set { platform_split_ch }

    TRIMMOMATIC(platform_split_ch.illumina)
    NANOFILT(platform_split_ch.ont)

    // --- 2. Alignment ---
    def trimmomatic_paired_reads = TRIMMOMATIC.out.reads
        .map { sample_id, platform, strain_id, r1p, r1u, r2p, r2u ->
            tuple(sample_id, platform, [r1p, r2p], strain_id)
        }

    def trimmed_reads_ch = trimmomatic_paired_reads.mix(NANOFILT.out.reads)
    
    def reads_for_alignment = trimmed_reads_ch
        .map { sample_id, platform, reads, strain_id ->
            def ref_fasta = file(params.references[strain_id])
            def ref_dir = ref_fasta.parent
            def ref_name = ref_fasta.name
            tuple(sample_id, platform, reads, ref_dir, ref_name, strain_id)
        }

    reads_for_alignment
        .branch { sample_id, platform, reads, ref_dir, ref_name, strain_id ->
            illumina: platform == 'illumina'
            ont:      platform == 'ont'
        }
        .set { align_split_ch }

    ALIGN_ILLUMINA(align_split_ch.illumina)
    ALIGN_ONT(align_split_ch.ont)

    // --- 3. Variant Calling ---
    def aligned_bams_ch = ALIGN_ILLUMINA.out.bam.mix(ALIGN_ONT.out.bam)

    def bams_for_variant_calling = aligned_bams_ch
        .map { sample_id, platform, bam, strain_id ->
            def ref_fasta = file(params.references[strain_id])
            tuple(sample_id, bam, ref_fasta, strain_id)
        }

    BCFTOOLS_MPILEUP(bams_for_variant_calling)
    BCFTOOLS_CONSENSUS(BCFTOOLS_MPILEUP.out.vcf)

    // --- 4. Nextclade Analysis ---
    EXTRACT_HA_SEGMENT(BCFTOOLS_CONSENSUS.out.fasta)
    
    def nextclade_input_ch = EXTRACT_HA_SEGMENT.out
        .map { sample_id, ha_fasta, strain_id ->
            def dataset_path = file(params.nextclade_datasets[strain_id])
            tuple(sample_id, ha_fasta, strain_id, dataset_path)
        }

    nextclade_results = NEXTCLADE(nextclade_input_ch)

def nextclade_csv = nextclade_results.csv
def nextclade_extra = nextclade_results.extra_outputs


    // --- 5. Final Reporting ---
    // Handle FASTQC output which is a list of files, and flatten it
    def fastqc_reports = FASTQC.out.reports.flatMap { strain_id, files ->
        files.collect { file -> tuple(strain_id, file) }
    }

    // All other report channels are already in the correct (strain_id, file) format
        def other_reports = TRIMMOMATIC.out.log
        .mix(
            ALIGN_ILLUMINA.out.log,
            ALIGN_ILLUMINA.out.flagstat,
            ALIGN_ONT.out.log,
            ALIGN_ONT.out.flagstat,
            BCFTOOLS_MPILEUP.out.stats,
            NEXTCLADE.out.csv
        )


    // Combine the FASTQC reports with the others
    def all_reports = fastqc_reports.mix(other_reports)

    // Group all the report files by strain_id. The output will be a channel of [strain_id, [file_list]]
    def multiqc_input_ch = all_reports.groupTuple()

    // Run MULTIQC once for each strain
    MULTIQC(multiqc_input_ch)
}