# scripts/extract_variants.R

library(vcfR)
library(dplyr)
library(tidyr)
library(stringr)

# Get list of all VCF files
vcf_files <- list.files("results/variants/h1n1/", pattern = "*.vcf.gz$", full.names = TRUE)

# Function to extract sample variants
extract_variants <- function(vcf_file) {
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  fix <- as.data.frame(vcf@fix)
  sample_name <- str_remove(basename(vcf_file), ".vcf.gz$")
  fix <- fix %>%
    mutate(Sample = sample_name,
           Pos = as.numeric(POS),
           Ref = REF,
           Alt = ALT,
           Variant_ID = paste0(CHROM, "_", POS, "_", REF, ">", ALT)) %>%
    select(Sample, Variant_ID, Pos, Ref, Alt)
  return(fix)
}

# Apply to all files and bind into one data frame
variant_list <- lapply(vcf_files, extract_variants)
all_variants <- bind_rows(variant_list)
