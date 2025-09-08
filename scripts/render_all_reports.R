library(rmarkdown)
library(dplyr)
library(here)

strains <- c("h1n1", "h3n2")

for (strain in strains) {
  message("Rendering overall report for ", strain)

  overall_report_dir <- here("results", "downstream", "reports", "overall")
  if (!dir.exists(overall_report_dir)) {
    dir.create(overall_report_dir, recursive = TRUE)
    if (!dir.exists(overall_report_dir)) {
      stop("Failed to create directory: ", overall_report_dir)
    }
  }

  rmd_file <- here("scripts", "strain_report.Rmd")

  # Render overall report
  render(
    input = rmd_file,
    params = list(strain_name = strain, sample_id = NULL),
    output_file = file.path(overall_report_dir, paste0(strain, "_overall_report.html")),
    knit_root_dir = here::here()  # ensures relative paths inside Rmd work
  )

  variants_path <- here("results", "downstream", "variants_bwa", strain, paste0(strain, "_variants.csv"))

  if (file.exists(variants_path)) {
    variants <- read.csv(variants_path)
    samples <- unique(variants$Sample)

    individual_report_dir <- here("results", "downstream", "reports", "individual")
    if (!dir.exists(individual_report_dir)) {
      dir.create(individual_report_dir, recursive = TRUE)
      if (!dir.exists(individual_report_dir)) {
        stop("Failed to create directory: ", individual_report_dir)
      }
    }

    for (sample_id in samples) {
      message("Rendering report for sample ", sample_id, " in strain ", strain)
      output_file <- file.path(individual_report_dir, paste0(sample_id, "_report.html"))

      render(
        input = rmd_file,
        params = list(strain_name = strain, sample_id = sample_id),
        output_file = output_file,
        knit_root_dir = here::here()
      )
    }
  } else {
    warning("No variant CSV found for strain ", strain)
  }
}
