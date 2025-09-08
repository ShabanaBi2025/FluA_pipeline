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
<<<<<<< HEAD

=======
  
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
  info_list <- lapply(info_split, function(x) {
    kv <- str_split(x, "=")
    sapply(kv, function(y) if(length(y) > 1) y[2] else NA)
  })
<<<<<<< HEAD

  info_df <- data.frame(matrix(NA, nrow = length(info_list), ncol = length(all_keys)))
  colnames(info_df) <- all_keys

=======
  
  info_df <- data.frame(matrix(NA, nrow = length(info_list), ncol = length(all_keys)))
  colnames(info_df) <- all_keys
  
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
  for (i in seq_along(info_list)) {
    vals <- info_list[[i]]
    names(vals) <- sapply(str_split(info_split[[i]], "="), `[`, 1)
    info_df[i, names(vals)] <- vals
  }
  info_df
}

# Extract variants and annotations from annotated VCF files
extract_variants_from_folder <- function(folder_path) {
<<<<<<< HEAD
  vcf_files <- list.files(folder_path, pattern = "\\.annotated\\.vcf$", full.names = TRUE, recursive = TRUE)

  print(paste("Found annotated VCF files:", paste(basename(vcf_files), collapse = ", ")))

  if (length(vcf_files) == 0) return(NULL)

  extract_variants <- function(vcf_file) {
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    fix <- as.data.frame(vcf@fix)
    sample_name <- str_remove(basename(vcf_file), "\\.annotated\\.vcf$")
    info_df <- extract_info(vcf)

    ann_df <- NULL
    if ("ANN" %in% colnames(info_df)) {
      ann_raw <- info_df$ANN
      ann_long <- unlist(str_split(ann_raw, ","))  # Flatten all annotations

      ann_per_variant <- sapply(str_split(ann_raw, ","), length)

      fix_expanded <- fix[rep(1:nrow(fix), times = ann_per_variant), ]
      sample_vec <- rep(sample_name, nrow(fix_expanded))

=======
  vcf_files <- list.files(folder_path, pattern = "\\.vcf.annotated\\.vcf$", full.names = TRUE, recursive = TRUE)
  if (length(vcf_files) == 0) return(NULL)
  
  extract_variants <- function(vcf_file) {
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    fix <- as.data.frame(vcf@fix)
    sample_name <- str_remove(basename(vcf_file), "\\.vcf.annotated\\.vcf$")
    info_df <- extract_info(vcf)
    
    ann_df <- NULL
    if ("ANN" %in% colnames(info_df)) {
      ann_raw <- info_df$ANN
      ann_long <- str_split(ann_raw, ",")
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
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
<<<<<<< HEAD

      ann_df$Sample <- sample_vec
      ann_df$Pos <- as.numeric(fix_expanded$POS)
      ann_df$CHROM <- fix_expanded$CHROM
      ann_df$Ref <- fix_expanded$REF
      ann_df$Alt <- fix_expanded$ALT

      ann_df$Variant_ID <- paste0(ann_df$CHROM, "_", ann_df$Pos, "_", ann_df$Ref, ">", ann_df$Alt)
    }

=======
      ann_df$Sample <- sample_name
    }
    
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
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
<<<<<<< HEAD

    list(variants = variants_basic, annotations = ann_df)
  }

  variant_data_list <- lapply(vcf_files, extract_variants)
  variants_all <- do.call(rbind, lapply(variant_data_list, `[[`, "variants"))
  annotations_all <- do.call(rbind, lapply(variant_data_list, `[[`, "annotations"))

=======
    
    list(variants = variants_basic, annotations = ann_df)
  }
  
  variant_data_list <- lapply(vcf_files, extract_variants)
  variants_all <- do.call(rbind, lapply(variant_data_list, `[[`, "variants"))
  annotations_all <- do.call(rbind, lapply(variant_data_list, `[[`, "annotations"))
  
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
  list(variants = variants_all, annotations = annotations_all)
}

# Summarize mutation annotations and save barplot
<<<<<<< HEAD
summarize_mutations <- function(annotations_df, output_file) {
=======
 
 summarize_mutations <- function(annotations_df, output_file) {
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
  if (is.null(annotations_df) || nrow(annotations_df) == 0) {
    message("No mutation annotations found")
    return(NULL)
  }
<<<<<<< HEAD

  # Split combined annotations separated by '&' into multiple rows
  annotations_df <- annotations_df %>%
    tidyr::separate_rows(Annotation, sep = "&")

  summary_table <- annotations_df %>%
    group_by(Annotation) %>%
    summarize(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  write.csv(summary_table, output_file, row.names = FALSE)
  message("Mutation summary saved to ", output_file)

  p <- ggplot(summary_table, aes(x = reorder(Annotation, -count), y = count, fill = Annotation)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.3, size = 3) +
    theme_classic() +
    labs(title = "Mutation Annotation Counts", x = "Annotation Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

=======
  
  summary_table <- annotations_df %>%
    group_by(Annotation) %>%
    summarize(count = n()) %>%
    arrange(desc(count))
  
  write.csv(summary_table, output_file, row.names = FALSE)
  message("Mutation summary saved to ", output_file)
  
  p <- ggplot(summary_table, aes(x = reorder(Annotation, -count), y = count, fill = Annotation)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.3, size = 3) +   # <-- Add count labels above bars
    theme_classic() +
    labs(title = "Mutation Annotation Counts", x = "Annotation Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
  ggsave(sub(".csv$", "_barplot.png", output_file), plot = p, width = 10, height = 6)
  message("Mutation summary bar plot saved to ", sub(".csv$", "_barplot.png", output_file))
}

<<<<<<< HEAD
# Plot variant distribution histogram by position
plot_variant_distribution <- function(variants_df, output_file) {
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    message("No variants data for distribution plot")
    return(NULL)
  }

  p <- ggplot(variants_df, aes(x = Pos)) +
    geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
    theme_classic() +
    labs(title = "Variant Position Distribution", x = "Genomic Position", y = "Count")

  ggsave(output_file, plot = p, width = 10, height = 6)
  message("Variant position distribution plot saved to ", output_file)
}

# Plot heatmap: Position vs Annotation Impact
plot_impact_heatmap <- function(annotations_df, output_file) {
  if (is.null(annotations_df) || nrow(annotations_df) == 0) {
    message("No annotations for impact heatmap")
    return(NULL)
  }

  annotations_df <- annotations_df %>%
    mutate(Pos_bin = floor(as.numeric(Pos) / 1000) * 1000)

  heatmap_df <- annotations_df %>%
    group_by(Pos_bin, Annotation_Impact) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(!is.na(Annotation_Impact))

  p <- ggplot(heatmap_df, aes(x = Annotation_Impact, y = Pos_bin, fill = count)) +
    geom_tile(color = "grey80") +
    scale_fill_gradient(low = "white", high = "red", na.value = "white") +
    theme_minimal() +
    labs(title = "Variant Heatmap (Position vs Annotation Impact)",
         x = "Annotation Impact", y = "Position Bin") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(output_file, plot = p, width = 10, height = 8)
  message("Position vs Annotation Impact heatmap saved to ", output_file)
}

# Pie chart of variant annotations (categories)
plot_annotation_pie <- function(annotations_df, output_file) {
  if (is.null(annotations_df) || nrow(annotations_df) == 0) {
    message("No annotations for annotation pie chart")
    return(NULL)
  }

  # Separate combined annotations joined by '&'
  annotations_df <- annotations_df %>%
    tidyr::separate_rows(Annotation, sep = "&")

  annotation_counts <- annotations_df %>%
    count(Annotation) %>%
    mutate(percentage = n / sum(n) * 100,
           label = paste0(Annotation, " (", round(percentage, 1), "%)"))

  p <- ggplot(annotation_counts, aes(x = "", y = n, fill = Annotation)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = "Variant Annotation Distribution") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )

  ggsave(output_file, plot = p, width = 6, height = 6)
  message("Annotation pie chart saved to ", output_file)
}

# Bar plot of annotation impact counts
plot_impact_barplot <- function(annotations_df, output_file) {
  if (is.null(annotations_df) || nrow(annotations_df) == 0) {
    message("No annotations for impact bar plot")
    return(NULL)
  }

  impact_summary <- annotations_df %>%
    count(Annotation_Impact) %>%
    filter(!is.na(Annotation_Impact))

  p <- ggplot(impact_summary, aes(x = reorder(Annotation_Impact, -n), y = n, fill = Annotation_Impact)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.3) +
    theme_classic() +
    labs(title = "Variant Impact Counts", x = "Annotation Impact", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(output_file, plot = p, width = 8, height = 6)
  message("Variant impact bar plot saved to ", output_file)
}

# NEW: Plot bar plot comparing variant annotation counts between strains
plot_annotation_comparison <- function(output_dir) {
  # SnpEff summary data for H1N1 and H3N2:
  variant_counts <- data.frame(
    Annotation = c("synonymous_variant", "missense_variant", "downstream_gene_variant",
                   "intron_variant", "intragenic_variant", "frameshift_variant", "splice_region_variant",
                   "intergenic_region", "upstream_gene_variant", "stop_gained", "stop_retained_variant", "stop_lost"),
    H1N1 = c(2682, 962, 412, 335, 85, 7, 5, 0, 0, 0, 0, 0),
    H3N2 = c(5761, 1846, 1976, 846, 0, 0, 21, 75, 49, 39, 7, 2)
  )

  variant_long <- variant_counts %>%
    tidyr::pivot_longer(cols = c(H1N1, H3N2), names_to = "Strain", values_to = "Count")

  p <- ggplot(variant_long, aes(x = Annotation, y = Count, fill = Strain)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_classic() +
    labs(title = "Variant Annotation Counts by Strain",
         x = "Annotation Type",
         y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_file <- file.path(output_dir, "variant_annotation_comparison.png")
  ggsave(output_file, plot = p, width = 12, height = 6)
  message("Variant annotation comparison plot saved to ", output_file)
=======

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
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
}

# Main processing function
process_all_strains <- function(input_dir = base_input_dir, output_dir = base_output_dir) {
  strain_folders <- list.dirs(input_dir, recursive = FALSE)
<<<<<<< HEAD

  if (length(strain_folders) == 0) {
    message("No strain folders found in ", input_dir)
    return(NULL)
  }

  for (strain_folder in strain_folders) {
    strain_name <- basename(strain_folder)
    message("Processing strain: ", strain_name)

=======
  
  for (strain_folder in strain_folders) {
    strain_name <- basename(strain_folder)
    message("Processing strain: ", strain_name)
    
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
    variant_data <- extract_variants_from_folder(strain_folder)
    if (is.null(variant_data)) {
      message("No annotated VCF files found for strain ", strain_name)
      next
    }
<<<<<<< HEAD

    variants_df <- variant_data$variants
    annotations_df <- variant_data$annotations

    if (!is.null(variants_df) && nrow(variants_df) > 0) {
      output_folder <- file.path(output_dir, strain_name)
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

      write.csv(variants_df, file.path(output_folder, paste0("variants_", strain_name, ".csv")), row.names = FALSE)
      message("Saved variants table for ", strain_name)

      if (!is.null(annotations_df) && nrow(annotations_df) > 0) {
        write.csv(annotations_df, file.path(output_folder, paste0("annotations_", strain_name, ".csv")), row.names = FALSE)
        message("Saved annotations table for ", strain_name)

        summarize_mutations(annotations_df, file.path(output_folder, paste0("mutation_summary_", strain_name, ".csv")))
        plot_variant_distribution(variants_df, file.path(output_folder, paste0("variant_distribution_", strain_name, ".png")))
        plot_impact_heatmap(annotations_df, file.path(output_folder, paste0("impact_heatmap_", strain_name, ".png")))
        plot_annotation_pie(annotations_df, file.path(output_folder, paste0("annotation_pie_", strain_name, ".png")))
        plot_impact_barplot(annotations_df, file.path(output_folder, paste0("impact_barplot_", strain_name, ".png")))
      }
    }
  }

  # Generate cross-strain variant annotation comparison plot
  plot_annotation_comparison(output_dir)
}

# Run the main process
process_all_strains()
this worked for that but then we changed it for the other plots fro   to this 
=======
    
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
>>>>>>> f502c5274a6e63da3e929c45feff7de4588a5fde
