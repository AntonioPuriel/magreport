# Regenerates the figures shown in README.md from the example output
# shipped with the package. Run with: source("data-raw/make_readme_figures.R")

devtools::load_all(".", quiet = TRUE)
library(ggplot2)

example_path <- function(file) {
  system.file("extdata", "example_run", file, package = "magreport")
}

mags    <- read_mag_quality(example_path("mag_quality.tsv"))
mapping <- read_mapping_summary(example_path("mapping_summary.tsv"))

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "man/figures/README-mag-quality.png",
  plot_mag_quality(mags, facet_by = "binner"),
  width = 8, height = 4, dpi = 150
)

ggsave(
  "man/figures/README-mapping-rates.png",
  plot_mapping_rates(mapping),
  width = 6, height = 2.6, dpi = 150
)

ggsave(
  "man/figures/README-taxonomy.png",
  plot_taxonomy(mags, rank = "genus", fill_by = "quality"),
  width = 6, height = 3, dpi = 150
)

message("Figures written to man/figures/")
