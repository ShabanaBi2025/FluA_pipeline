// main.nf
nextflow.enable.dsl = 2

// Only include necessary modules for QC and Trimming
include { fastqc }      from './modules/qc.nf'
include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { NANOFILT }    from './modules/nanofilt.nf'
// Remove includes for BWA_INDEX, ALIGN_ILLUMINA, ALIGN_ONT, BCFTOOLS_CALL

workflow {

  // --- 1. Create Raw Reads Channel - Grouping by Sample ID and Type ---
  Channel
    .fromPath(params.raw_reads) // Assumes params.raw_reads is set in base.config
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

  // Debugging: What does all_raw_samples emit before FastQC?
  all_raw_samples.subscribe { item -> println "DEBUG_ALL_RAW_SAMPLES: $item" }

  // --- Part 1: QC and Trimming ---
  fastqc_output = fastqc(all_raw_samples)

  qc_branched_output = fastqc_output
    .branch { sample_id, platform, reads, strain_id ->
      illumina: platform == 'illumina'
      ont:      platform == 'ont'
    }

  // Debugging: What does each branch emit for trimming?
  qc_branched_output.illumina.subscribe { item -> println "DEBUG_TRIM_INPUT: Illumina branch for trimming: $item" }
  qc_branched_output.ont.subscribe { item -> println "DEBUG_TRIM_INPUT: ONT branch for trimming: $item" }

  // TRIMMOMATIC for Illumina data
  // TRIMMOMATIC outputs a flat tuple: (sample_id, platform, strain_id, paired1_str, unpaired1_str, paired2_str, unpaired2_str)
  def trimmomatic_output_channels = TRIMMOMATIC(qc_branched_output.illumina)
  
  // Map Trimmomatic's output to the desired (sample_id, platform, [paired_reads_list], strain_id) format
  // This is 'trimmed_illumina_reads'
  trimmed_illumina_reads = trimmomatic_output_channels
    .map { sample_id, platform, strain_id, paired1_str, unpaired1_str, paired2_str, unpaired2_str ->
        def paired_reads_list = [file(paired1_str), file(paired2_str)] // Convert string paths to Path objects
        tuple(sample_id, platform, paired_reads_list, strain_id)
    }
  // Debugging: What does trimmed_illumina_reads emit?
  trimmed_illumina_reads.subscribe { item -> println "DEBUG_TRIMMED_ILLUMINA_OUTPUT: $item" }


  // NANOFILT for ONT data
  // NANOFILT outputs: (sample_id, platform, nanofilt_read_path, strain_id)
  filtered_ont_reads = NANOFILT(qc_branched_output.ont)
  // Debugging: What does filtered_ont_reads emit?
  filtered_ont_reads.subscribe { item -> println "DEBUG_NANOFILT_OUTPUT: $item" }

  // --- End of Part 1 (QC and Trimming) ---
  // No further steps in this simplified workflow.
}