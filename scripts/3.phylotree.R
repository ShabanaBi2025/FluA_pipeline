library(ape)
library(ggtree)
library(dplyr)
library(readr)
library(ggplot2)

# List of strains you want to process
strains <- c("h1n1", "h3n2")

for (strain in strains) {
  message("Processing strain: ", strain)

  # Define input paths
  base_dir <- file.path("results", "nextclade", strain, paste0(strain, "_nextclade_results"))
  tree_file <- file.path(base_dir, "nextclade.nwk")
  metadata_file <- file.path(base_dir, "nextclade.tsv")

  # Validate files exist
  if (!file.exists(tree_file) || !file.exists(metadata_file)) {
    warning("Missing files for strain: ", strain)
    next
  }

  # Fix tree if missing semicolon
  tree_text <- readLines(tree_file)
  if (!grepl(";$", tree_text[length(tree_text)])) {
    cat("Appending semicolon to tree file for", strain, "\n")
    writeLines(c(tree_text, ";"), tree_file)
  }

  # Load tree and metadata
  tree <- read.tree(tree_file)
  metadata <- read_tsv(metadata_file, show_col_types = FALSE)

  # Rename seqName to sample_id if needed
  if ("seqName" %in% names(metadata)) {
    metadata <- metadata %>% rename(sample_id = seqName)
  }

  # Join tree tip labels with metadata
  tree_data <- full_join(
    tibble(sample_id = tree$tip.label),
    metadata,
    by = "sample_id"
  )

  # Plot tree colored by clade
  p <- ggtree(tree, size = 0.25, color = "grey30") %<+% tree_data +
    geom_tiplab(aes(color = as.factor(clade)), size = 3) +
    scale_color_brewer(palette = "Set3") +
    theme_tree2() +
    labs(title = paste0("Phylogenetic Tree - ", toupper(strain))) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10)
    )

  # Output directory and save plot
  output_dir <- file.path("results", "downstream", "phylogenetic", strain)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(file.path(output_dir, paste0(strain, "_tree_plot.png")), plot = p, width = 10, height = 8)

  message("Saved tree plot for: ", strain, "\n")
}
