// main.nf
nextflow.enable.dsl = 2

include { fastqc }      from './modules/qc.nf'
include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { NANOFILT }    from './modules/nanofilt.nf'
include { BWA_INDEX }   from './modules/bwa_index.nf'
include { ALIGN_ILLUMINA } from './modules/align_illumina.nf'
include { ALIGN_ONT }      from './modules/align_ont.nf'

workflow {

  // --- 1. Define Reference Genome Channel and Index ---
  Channel
    .from(params.references.entrySet())
    .map { entry -> tuple(entry.key, file(entry.value)) }
    .set { all_reference_fastas }

  indexed_references_channel = BWA_INDEX(all_reference_fastas)

  // DEBUG: What does indexed_references_channel emit? (Key for join is element 0)
  indexed_references_channel.subscribe { item -> println "DEBUG_JOIN_INPUT: indexed_references_channel emitted: $item" }


  // --- 2. Create Raw Reads Channel - Grouping by Sample ID and Type ---
  Channel
    .fromPath(params.raw_reads)
    .map { file ->
            def strain_id = file.parent.name
            def filename = file.baseName
            def sample_id_clean = filename.replaceAll(/_R?[12]|\.fastq(.gz)?/, '')
            def platform_type
            if (filename.endsWith('_1') || filename.endsWith('_R1') ||
                filename.endsWith('_2') || filename.endsWith('_R2')) {
                platform_type = 'illumina'
            } else {
                platform_type = 'ont'
            }
            return tuple(strain_id, sample_id_clean, platform_type, file)
        }
        .groupTuple(by: [0, 1])
        .map { it ->
            def current_strain_id = it.get(0)
            def current_sample_id_clean = it.get(1)
            def platform_types_list = it.get(2)
            def files_list = it.get(3)

            def final_platform_type = platform_types_list.get(0)
            def final_files = files_list
            if (final_platform_type == 'illumina') {
                final_files.sort { a, b -> a.name.compareTo(b.name) }
            }
            return tuple(current_sample_id_clean, final_platform_type, final_files, current_strain_id)
        }
        .set { all_raw_samples }

  fastqc_output = fastqc(all_raw_samples)

  qc_branched_output = fastqc_output
    .branch { sample_id, platform, reads, strain_id ->
      illumina: platform == 'illumina'
      ont:      platform == 'ont'
    }

  def trimmomatic_output_channels = TRIMMOMATIC(qc_branched_output.illumina)
  trimmed_metadata_channel = trimmomatic_output_channels[0]
  trimmed_files_channel = trimmomatic_output_channels[1]

  trimmed_illumina_reads = trimmed_metadata_channel
    .join(trimmed_files_channel, by: 0)
    .map { sample_id, platform, strain_id, paired1, unpaired1, paired2, unpaired2 ->
        tuple(sample_id, platform, [paired1, paired2], strain_id)
    }

  filtered_ont_reads = NANOFILT(qc_branched_output.ont)

  channel_for_align_illumina = trimmed_illumina_reads
  
  channel_for_align_ont = filtered_ont_reads
    .map { sample_id, platform, single_read, strain_id ->
      tuple(sample_id, platform, [single_read], strain_id)
    }

  final_reads_for_align = channel_for_align_illumina.mix(channel_for_align_ont)

  // DEBUG: What does final_reads_for_align emit? (Key for join is element 3)
  final_reads_for_align.subscribe { item -> println "DEBUG_JOIN_INPUT: final_reads_for_align emitted: $item" }


  // --- 3. JOIN Reads with Indexed References for Alignment ---
  align_inputs = final_reads_for_align
    .join(indexed_references_channel, by: 3) // Join by strain_id from reads (index 3), and ref_id from indexed_references (index 0)
    .map { sample_id, platform, reads_list, read_strain_id, indexed_ref_id, fasta_file, amb, ann, bwt, pac, sa ->
      // Pass the strain_id through to ALIGN's input, so ALIGN can output it
      tuple(sample_id, platform, reads_list, fasta_file, read_strain_id)
    }
    .set { align_inputs_with_ref }


  // --- 4. Branch Aligned Inputs by Platform and Send to Specific Aligners ---
  align_inputs_with_ref
    .branch { sample_id, platform, reads_list, fasta_file, strain_id ->
      illumina_align: platform == 'illumina'
      ont_align:      platform == 'ont'
    }
    .set { align_branched_output }

  // Call specific alignment processes (ONLY ONCE)
  ALIGN_ILLUMINA(align_branched_output.illumina_align)
  ALIGN_ONT(align_branched_output.ont_align)

  // Output for downstream (e.g., variant calling) would then combine ALIGN_ILLUMINA.out and ALIGN_ONT.out
  // E.g.: final_aligned_bams = ALIGN_ILLUMINA.out.mix(ALIGN_ONT.out)
  // BCFTOOLS_CALL(final_aligned_bams)
}