// test_input_channel.nf
nextflow.enable.dsl = 2

// Parameters are defined as simple Groovy variables
def TRIMMED_ILLUMINA_READS_PATTERN_GENERIC = "results/trimmomatic/**/*.fq.gz"
def NANOFILT_ONT_READS_PATTERN = "results/nanofilt/**/*.fastq.gz" // ONT part is working, keep as is

workflow {
  // --- 1. Test Indexed References Channel (Unchanged, already working) ---
  def indexed_references_channel = Channel.from(
    tuple('h1n1',
          file('results/indexed_reference/h1n1/h1n1_reference.fasta'),
          file('results/indexed_reference/h1n1/h1n1_reference.fasta.amb'),
          file('results/indexed_reference/h1n1/h1n1_reference.fasta.ann'),
          file('results/indexed_reference/h1n1/h1n1_reference.fasta.bwt'),
          file('results/indexed_reference/h1n1/h1n1_reference.fasta.pac'),
          file('results/indexed_reference/h1n1/h1n1_reference.fasta.sa')),
    tuple('h3n2',
          file('results/indexed_reference/h3n2/h3n2_reference.fasta'),
          file('results/indexed_reference/h3n2/h3n2_reference.fasta.amb'),
          file('results/indexed_reference/h3n2/h3n2_reference.fasta.ann'),
          file('results/indexed_reference/h3n2/h3n2_reference.fasta.bwt'),
          file('results/indexed_reference/h3n2/h3n2_reference.fasta.pac'),
          file('results/indexed_reference/h3n2/h3n2_reference.fasta.sa'))
  )
  indexed_references_channel.subscribe { item -> println "TEST_OUTPUT_REFERENCES: Emitted item: $item" }


  // --- 2. Test Input Reads Channels (Dynamic Pairing for Illumina with fromPath) ---
  // Channel for ALIGN_ILLUMINA (from Trimmomatic output)
  Channel
    .fromPath(TRIMMED_ILLUMINA_READS_PATTERN_GENERIC)
    .map { file ->
        def full_path = file.toString()
        def filename = file.name // e.g., ERR4552844_1_paired.fq.gz

        println "DEBUG_ILLUMINA_PARSE: Processing file: $file"
        println "DEBUG_ILLUMINA_PARSE: file.name (raw): '$filename'"
        println "DEBUG_ILLUMINA_PARSE: file.name (bytes): '${filename.bytes.encodeHex()}'" // Crucial diagnostic

        def strain_id
        if (full_path.contains('/h1n1/')) { strain_id = 'h1n1' }
        else if (full_path.contains('/h3n2/')) { strain_id = 'h3n2' }
        else { strain_id = 'unknown_strain' }

        def sample_id_base = null
        def read_category = null
        
        // --- CRITICAL CHANGE: Primitive extraction based on fixed suffixes ---
        // Try to replace the known suffixes to get the sample ID base
        // Then, determine category from what's left of the suffix.
        if (filename.endsWith('.fq.gz')) { // Check for the overall extension
            def base_name_without_ext = filename.substring(0, filename.length() - '.fq.gz'.length()) // e.g. ERR4552844_1_paired

            if (base_name_without_ext.endsWith('_1_paired')) {
                sample_id_base = base_name_without_ext.replace('_1_paired', '')
                read_category = 'paired_r1'
            } else if (base_name_without_ext.endsWith('_2_paired')) {
                sample_id_base = base_name_without_ext.replace('_2_paired', '')
                read_category = 'paired_r2'
            } else if (base_name_without_ext.endsWith('_1_unpaired') || base_name_without_ext.endsWith('_2_unpaired')) {
                // This correctly filters out unpaired reads
                println "WARN: Found unpaired Illumina read (matched by primitive suffix check): $file. This will be skipped for paired alignment."
                return null
            } else {
                // This 'else' means it's a .fq.gz file but doesn't have a recognizable _paired or _unpaired suffix.
                println "WARN: Skipping non-standard Illumina file (no recognizable _paired/_unpaired suffix after primitive check): $file"
                return null
            }
        } else {
            // This 'else' block catches files not ending in .fq.gz (e.g., .fastq.gz for ONT, or others)
            // It assumes ONT files will be handled by the ONT channel (which is working).
            // For this Illumina channel, it filters out non-.fq.gz files.
            println "WARN: Skipping non-standard Illumina file (not ending in .fq.gz): $file"
            return null
        }

        // Output tuple for grouping
        tuple(strain_id, sample_id_base, read_category, file)
    }
    .filter { it != null } // Filter out nulls from the map operation (non-standard files)
    .groupTuple(by: [0, 1]) // Group by strain_id and sample_id_base
    .map { grouped_keys, files_with_categories ->
        def current_strain_id = grouped_keys.get(0)
        def current_sample_id_base = grouped_keys.get(1)

        def illumina_r1 = null
        def illumina_r2 = null

        files_with_categories.each { item_tuple ->
            def category = item_tuple.get(0) // This is "paired_r1" or "paired_r2"
            def file_path = item_tuple.get(1)
            if (category == 'paired_r1') {
                illumina_r1 = file_path
            } else if (category == 'paired_r2') {
                illumina_r2 = file_path
            }
        }

        if (illumina_r1 && illumina_r2) {
            def final_files = [illumina_r1, illumina_r2].sort { a, b -> a.name.compareTo(b.name) }
            return tuple(current_sample_id_base, 'illumina', final_files, current_strain_id)
        } else {
            println "WARN: Found unmatched Illumina pair for sample ${current_sample_id_base} in strain ${current_strain_id}. Skipping for alignment."
            return null
        }
    }
    .filter { it != null } // Filter out any nulls from the final map (unmatched pairs)
    .set { channel_for_align_illumina }
    
  channel_for_align_illumina.subscribe { item -> println "TEST_OUTPUT_ILLUMINA: Emitted item: $item" }


  // Channel for ALIGN_ONT (from NanoFilt output) - Unchanged, already working
  Channel
    .fromPath(NANOFILT_ONT_READS_PATTERN)
    .map { file ->
        def sample_id_clean = file.baseName.replaceAll(/_nanofilt|\.fastq\.gz/, '')
        def strain_id
        def full_path = file.toString()
        if (full_path.contains('/h1n1/')) {
            strain_id = 'h1n1'
        } else if (full_path.contains('/h3n2/')) {
            strain_id = 'h3n2'
        } else {
            strain_id = 'unknown_strain'
        }
        tuple(sample_id_clean, 'ont', [file], strain_id)
    }
    .subscribe { item -> println "TEST_OUTPUT_ONT: Emitted item: $item" }

  // No processes will be run, only channel operations and debug printing
}