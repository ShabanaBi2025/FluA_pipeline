// test_string_compare.nf
nextflow.enable.dsl = 2

// Hardcode one problematic filename for precise testing
def TEST_ILLUMINA_FILE_PATH = "results/trimmomatic/h1n1/illumina/ERR4552844_1_paired.fq.gz"
def SUFFIX_TO_CHECK = "_1_paired.fq.gz"

workflow {
  Channel
    .fromPath(TEST_ILLUMINA_FILE_PATH)
    .map { file ->
        def filename = file.name
        def contains_result = filename.contains(SUFFIX_TO_CHECK)

        println "DIAGNOSTIC: --- File String Check ---"
        println "DIAGNOSTIC: Full path: $file"
        println "DIAGNOSTIC: Filename as read (file.name): '$filename'"
        println "DIAGNOSTIC: Suffix to check: '$SUFFIX_TO_CHECK'"
        println "DIAGNOSTIC: filename.contains(suffix_to_check) result: $contains_result"

        println "DIAGNOSTIC: --- Character-by-Character Comparison ---"
        println "DIAGNOSTIC: Filename length: ${filename.length()}"
        println "DIAGNOSTIC: Suffix length: ${SUFFIX_TO_CHECK.length()}"

        // Print character by character with ASCII/Unicode values
        for (i = 0; i < filename.length(); i++) {
            def char_file = filename.charAt(i)
            def ascii_file = (int)char_file
            def hex_file = String.format("%02X", ascii_file)
            println "DIAGNOSTIC: Filename char[$i]: '$char_file' (ASCII: $ascii_file, Hex: $hex_file)"
        }
        for (i = 0; i < SUFFIX_TO_CHECK.length(); i++) {
            def char_suffix = SUFFIX_TO_CHECK.charAt(i)
            def ascii_suffix = (int)char_suffix
            def hex_suffix = String.format("%02X", ascii_suffix)
            println "DIAGNOSTIC: Suffix char[$i]: '$char_suffix' (ASCII: $ascii_suffix, Hex: $hex_suffix)"
        }

        if (!contains_result) {
            println "DIAGNOSTIC: !!! CONTAINS FAILED - MANUAL CHECK NEEDED !!!"
        }
        
        return "Done" // Just return something to keep channel active
    }
    .subscribe { it -> } // Consume the channel
}