library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# Base input and output directories
base_input_dir <- "results/variants/"
base_output_dir <- "results/downstream/"

# Function to extract and classify variants from all VCFs in a given folder
# FILTERING to keep only SNPs
extract_variants_from_folder <- function(folder_path) {
  vcf_files <- list.files(folder_path, pattern = "\\.vcf.gz$", full.names = TRUE)
  if (length(vcf_files) == 0) return(NULL)  # No VCFs found
  
  extract_variants <- function(vcf_file) {
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    fix <- as.data.frame(vcf@fix)
    sample_name <- str_remove(basename(vcf_file), ".vcf.gz$")
    fix %>%
      mutate(
        Sample = sample_name,
        Pos = as.numeric(POS),
        Ref = REF,
        Alt = ALT,
        Variant_ID = paste0(CHROM, "_", POS, "_", REF, ">", ALT),
        Type = case_when(
          nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNP",
          nchar(REF) > nchar(ALT) ~ "Deletion",
          nchar(REF) < nchar(ALT) ~ "Insertion",
          TRUE ~ "Other"
        )
      ) %>%
      filter(Type == "SNP") %>%   # KEEP ONLY SNPs
      select(Sample, Variant_ID, Pos, Ref, Alt, Type)
  }
  
  variant_list <- lapply(vcf_files, extract_variants)
  bind_rows(variant_list)
}

# Plot variant distribution (bar plot)
plot_variant_distribution <- function(variants, output_file) {
  p <- variants %>%
    count(Pos) %>%
    ggplot(aes(x = Pos, y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Variant Counts Across Genome Positions", x = "Position", y = "Variant Count") +
    theme_classic()
  
  ggsave(output_file, plot = p, width = 8, height = 4)
}

# Plot heatmap: variant presence/absence across samples and positions
plot_variant_heatmap <- function(variants, output_file) {
  variant_matrix <- variants %>%
    mutate(Present = 1L) %>%
    select(Sample, Pos, Present) %>%
    distinct(Sample, Pos, .keep_all = TRUE)
  
  variant_wide <- pivot_wider(
    variant_matrix,
    names_from = Sample,
    values_from = Present,
    values_fill = 0L
  )
  
  variant_long <- pivot_longer(variant_wide, -Pos, names_to = "Sample", values_to = "Present")
  
  p <- ggplot(variant_long, aes(x = Sample, y = factor(Pos), fill = factor(Present))) +
    geom_tile(color = "grey70") +
    scale_fill_manual(values = c("0" = "white", "1" = "firebrick"), name = "Variant Present") +
    labs(title = "Variant Presence/Absence Heatmap",
         x = "Sample",
         y = "Genome Position") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(output_file, plot = p, width = 10, height = 8)
}

# Modified pie chart: adapts if only SNPs present
plot_variant_type_pie <- function(variants, output_file) {
  type_counts <- variants %>%
    count(Type) %>%
    mutate(prop = n / sum(n) * 100)
  
  if (nrow(type_counts) == 1 && type_counts$Type == "SNP") {
    p <- ggplot(type_counts, aes(x = Type, y = prop, fill = Type)) +
      geom_col() +
      geom_text(aes(label = paste0(round(prop), "%")), vjust = -0.5) +
      labs(title = "Variant Type Composition", y = "Percentage", x = NULL) +
      theme_minimal() +
      theme(legend.position = "none")
  } else {
    p <- ggplot(type_counts, aes(x = "", y = prop, fill = Type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      geom_text(aes(label = paste0(round(prop), "%")),
                position = position_stack(vjust = 0.5)) +
      labs(title = "Proportion of Variant Types", fill = "Variant Type", x = NULL, y = NULL) +
      theme_void()
  }
  
  ggsave(output_file, plot = p, width = 6, height = 5)
}

# Main processing function
process_all_strains <- function(input_dir = base_input_dir, output_dir = base_output_dir) {
  strain_folders <- list.dirs(input_dir, recursive = FALSE)
  
  for (strain_folder in strain_folders) {
    strain_name <- basename(strain_folder)
    message("Processing strain: ", strain_name)
    
    variants_df <- extract_variants_from_folder(strain_folder)
    
    if (!is.null(variants_df) && nrow(variants_df) > 0) {
      output_folder <- file.path(output_dir, strain_name)
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      
      # Save CSV
      output_file <- file.path(output_folder, paste0(strain_name, "_variants.csv"))
      write.csv(variants_df, output_file, row.names = FALSE)
      message("Saved variants for ", strain_name, " to ", output_file)
      
      # Generate plots
      plot_variant_distribution(variants_df, file.path(output_folder, "variant_distribution.png"))
      plot_variant_heatmap(variants_df, file.path(output_folder, "variant_heatmap.png"))
      
      # Only plot pie chart if more than one variant type present
      variant_types <- unique(variants_df$Type)
      if (length(variant_types) > 1) {
        plot_variant_type_pie(variants_df, file.path(output_folder, "variant_type_pie.png"))
        message("Pie chart saved for ", strain_name)
      } else {
        message("Pie chart skipped for ", strain_name, " (only one variant type present)")
      }
      
      message("Plots saved for ", strain_name)
    } else {
      message("No SNP variants found for strain ", strain_name)
    }
  }
}

# Run the script
process_all_strains()
