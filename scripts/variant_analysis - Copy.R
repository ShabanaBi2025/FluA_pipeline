library(vcfR)
library(dplyr)
library(stringr)
library(ggplot2)

# Base input and output directories
base_input_dir <- "results/variants"
base_output_dir <- "results/downstream"

# Function to extract variants from all VCFs in a given folder
extract_variants_from_folder <- function(folder_path) {
  vcf_files <- list.files(folder_path, pattern = "\\.vcf\\.gz$", full.names = TRUE)
  if (length(vcf_files) == 0) return(NULL)  # No VCF files found
  
  extract_variants <- function(vcf_file) {
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    fix <- as.data.frame(vcf@fix)
    sample_name <- str_remove(basename(vcf_file), "\\.vcf\\.gz$")
    fix %>%
      mutate(Sample = sample_name,
             Pos = as.numeric(POS),
             Ref = REF,
             Alt = ALT,
             Variant_ID = paste0(CHROM, "_", POS, "_", REF, ">", ALT)) %>%
      select(Sample, Variant_ID, Pos, Ref, Alt)
  }
  
  variant_list <- lapply(vcf_files, extract_variants)
  bind_rows(variant_list)
}

# Function to plot variant distribution
plot_variant_distribution <- function(variants, output_file) {
  p <- variants %>%
    count(Pos) %>%
    ggplot(aes(x = Pos, y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Variant Counts Across Genome Positions", x = "Position", y = "Variant Count") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black", size = 14, face = "bold")
    )
  
  ggsave(output_file, plot = p, width = 8, height = 4, bg = "white")
}


# Main loop to process each strain folder
strain_folders <- list.dirs(base_input_dir, recursive = FALSE)

for (strain_folder in strain_folders) {
  strain_name <- basename(strain_folder)
  message("Processing strain: ", strain_name)
  
  variants_df <- extract_variants_from_folder(strain_folder)
  
  if (!is.null(variants_df)) {
    # Create output folder if needed
    output_folder <- file.path(base_output_dir, strain_name)
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
    
    # Save CSV file
    output_csv <- file.path(output_folder, paste0(strain_name, "_variants.csv"))
    write.csv(variants_df, output_csv, row.names = FALSE)
    message("Saved variants CSV: ", output_csv)
    
    # Save plot
    output_plot <- file.path(output_folder, paste0(strain_name, "_variant_distribution.png"))
    plot_variant_distribution(variants_df, output_plot)
    message("Saved variant distribution plot: ", output_plot)
  } else {
    message("No VCF files found for strain: ", strain_name)
  }
}
