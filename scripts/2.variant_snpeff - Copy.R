library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# Base input and output directories
base_input_dir <- "results/snpeff/output/"
base_output_dir <- "results/downstream/variant_snpeff/"

# Function to extract INFO fields as dataframe from VCF object
extract_info <- function(vcf) {
  info_raw <- vcf@fix[, "INFO"]
  info_split <- str_split(info_raw, ";")
  all_keys <- unique(unlist(lapply(info_split, function(x) str_split(x, "=") %>% sapply(`[`, 1))))
  
  info_list <- lapply(info_split, function(x) {
    kv <- str_split(x, "=")
    sapply(kv, function(y) if(length(y) > 1) y[2] else NA)
  })
  
  info_df <- data.frame(matrix(NA, nrow = length(info_list), ncol = length(all_keys)))
  colnames(info_df) <- all_keys
  
  for (i in seq_along(info_list)) {
    vals <- info_list[[i]]
    names(vals) <- sapply(str_split(info_split[[i]], "="), `[`, 1)
    info_df[i, names(vals)] <- vals
  }
  info_df
}

# Extract variants and annotations from annotated VCF files
extract_variants_from_folder <- function(folder_path) {
  # Fixed pattern here to match *.annotated.vcf files
  vcf_files <- list.files(folder_path, pattern = "\\.annotated\\.vcf$", full.names = TRUE, recursive = TRUE)
  
  print(paste("Found annotated VCF files:", paste(basename(vcf_files), collapse = ", ")))  # Debug print
  
  if (length(vcf_files) == 0) return(NULL)
  
  extract_variants <- function(vcf_file) {
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    fix <- as.data.frame(vcf@fix)
    sample_name <- str_remove(basename(vcf_file), "\\.annotated\\.vcf$")
    info_df <- extract_info(vcf)
    
    ann_df <- NULL
    if ("ANN" %in% colnames(info_df)) {
      ann_raw <- info_df$ANN
      ann_long <- str_split(ann_raw, ",")
      ann_list <- lapply(ann_long, function(x) {
        fields <- str_split(x, "\\|")[[1]]
        length(fields) <- 16
        fields
      })
      ann_df <- do.call(rbind, ann_list) %>% as.data.frame(stringsAsFactors = FALSE)
      colnames(ann_df) <- c("Allele","Annotation","Annotation_Impact","Gene_Name",
                            "Gene_ID","Feature_Type","Feature_ID","Transcript_BioType",
                            "Rank","HGVS.c","HGVS.p","cDNA.pos / cDNA.length",
                            "CDS.pos / CDS.length","AA.pos / AA.length","Distance","Errors")
      ann_df$Sample <- sample_name
    }
    
    variants_basic <- fix %>%
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
      )
    
    list(variants = variants_basic, annotations = ann_df)
  }
  
  variant_data_list <- lapply(vcf_files, extract_variants)
  variants_all <- do.call(rbind, lapply(variant_data_list, `[[`, "variants"))
  annotations_all <- do.call(rbind, lapply(variant_data_list, `[[`, "annotations"))
  
  list(variants = variants_all, annotations = annotations_all)
}

# Summarize mutation annotations and save barplot
summarize_mutations <- function(annotations_df, output_file) {
  if (is.null(annotations_df) || nrow(annotations_df) == 0) {
    message("No mutation annotations found")
    return(NULL)
  }
  
  summary_table <- annotations_df %>%
    group_by(Annotation) %>%
    summarize(count = n()) %>%
    arrange(desc(count))
  
  write.csv(summary_table, output_file, row.names = FALSE)
  message("Mutation summary saved to ", output_file)
  
  p <- ggplot(summary_table, aes(x = reorder(Annotation, -count), y = count, fill = Annotation)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.3, size = 3) +
    theme_classic() +
    labs(title = "Mutation Annotation Counts", x = "Annotation Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(sub(".csv$", "_barplot.png", output_file), plot = p, width = 10, height = 6)
  message("Mutation summary bar plot saved to ", sub(".csv$", "_barplot.png", output_file))
}

# Plot distribution of variants by genome position
plot_variant_distribution <- function(variants_df, output_file) {
  p <- ggplot(variants_df, aes(x = Pos)) +
    geom_histogram(binwidth = 100, fill = "steelblue", color = "black") +
    theme_classic() +
    labs(title = "Variant Position Distribution", x = "Position", y = "Count")
  
  ggsave(output_file, plot = p, width = 10, height = 6)
  message("Variant distribution plot saved to ", output_file)
}

# Plot heatmap: Position vs Variant Type
plot_variant_type_heatmap <- function(variants_df, output_file) {
  heatmap_df <- variants_df %>%
    count(Pos, Type) %>%
    complete(Pos = full_seq(Pos, 1), Type, fill = list(n = 0))
  
  p <- ggplot(heatmap_df, aes(x = Type, y = Pos, fill = n)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    labs(title = "Variant Heatmap (Position vs Type)", x = "Variant Type", y = "Position")
  
  ggsave(output_file, plot = p, width = 8, height = 10)
  message("Position vs Type heatmap saved to ", output_file)
}

# Pie chart of variant types
plot_variant_type_pie <- function(variants_df, output_file) {
  variant_counts <- variants_df %>%
    count(Type) %>%
    mutate(percentage = n / sum(n) * 100,
           label = paste0(Type, " (", round(percentage, 1), "%)"))
  
  p <- ggplot(variant_counts, aes(x = "", y = n, fill = Type)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = "Variant Type Distribution")
  
  ggsave(output_file, plot = p, width = 6, height = 6)
  message("Variant type pie chart saved to ", output_file)
}

# Main processing function
process_all_strains <- function(input_dir = base_input_dir, output_dir = base_output_dir) {
  strain_folders <- list.dirs(input_dir, recursive = FALSE)
  
  for (strain_folder in strain_folders) {
    strain_name <- basename(strain_folder)
    message("Processing strain: ", strain_name)
    
    variant_data <- extract_variants_from_folder(strain_folder)
    if (is.null(variant_data)) {
      message("No annotated VCF files found for strain ", strain_name)
      next
    }
    
    variants_df <- variant_data$variants
    annotations_df <- variant_data$annotations
    
    if (!is.null(variants_df) && nrow(variants_df) > 0) {
      output_folder <- file.path(output_dir, strain_name)
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      
      # Save variants CSV
      variants_file <- file.path(output_folder, paste0(strain_name, "_variants.csv"))
      write.csv(variants_df, variants_file, row.names = FALSE)
      message("Saved variants for ", strain_name, " to ", variants_file)
      
      # Summarize mutation annotations and plot barplot
      summarize_mutations(annotations_df, file.path(output_folder, paste0(strain_name, "_mutation_summary.csv")))
      
      # Plot variant position distribution
      plot_variant_distribution(variants_df, file.path(output_folder, "variant_distribution.png"))
      
      # Plot variant type heatmap (position vs type)
      plot_variant_type_heatmap(variants_df, file.path(output_folder, "variant_type_heatmap.png"))
      
      # Pie chart of variant types (only if more than one type)
      if (length(unique(variants_df$Type)) > 1) {
        plot_variant_type_pie(variants_df, file.path(output_folder, "variant_type_pie.png"))
        message("Pie chart saved for ", strain_name)
      } else {
        message("Pie chart skipped for ", strain_name, " (only one variant type present)")
      }
      
      message("Plots saved for ", strain_name)
    } else {
      message("No variants found for strain ", strain_name)
    }
  }
}

# Run the script
process_all_strains()
