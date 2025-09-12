library(ape)
library(ggtree)
library(dplyr)
library(readr)
library(ggplot2)

# List of strains to process
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

  #  tree 
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

  # Create a numeric label for each tip for clean plotting
  label_map <- tibble(
    sample_id = tree$tip.label,
    label_num = seq_along(tree$tip.label)
  )

  # Save label map as CSV for legend/reference
  output_dir <- file.path("results", "downstream", "phylogenetic", strain)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(label_map, file.path(output_dir, paste0(strain, "_label_legend.csv")))

  # Join metadata with label map
  metadata <- metadata %>% left_join(label_map, by = "sample_id")

  # Join tree tips with metadata including label_num
  tree_data <- full_join(
    tibble(sample_id = tree$tip.label),
    metadata,
    by = "sample_id"
  )

  # Plot tree with numeric labels and clade colors
 p <- ggtree(tree, size = 0.5, color = "grey80") %<+% tree_data +
  geom_tiplab(
    aes(label = label_num, color = as.factor(clade)),
    size = 3.5,
    align = FALSE,    # flush labels
    offset = 0        # no extra space, no leader lines
  ) +
scale_color_manual(
  values = c(
    "#E41A1C",  # Red
    "#377EB8",  # Blue
    "#984EA3",  # Purple
    "#4DAF4A",  # Green
    "#FF7F00",  # Orange
    "#A65628",  # Brown
    "#F781BF",  # Pink
    "#999999"   # Grey
  )
) +
  theme_tree2() +
  labs(title = paste0("Phylogenetic Tree - ", toupper(strain))) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

  # Save the plot
  ggsave(file.path(output_dir, paste0(strain, "_tree_plot.png")), plot = p, width = 10, height = 8, dpi = 300)

  message("Saved tree plot and label legend for: ", strain, "\n")
}
