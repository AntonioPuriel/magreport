example_path <- function(file) {
  system.file("extdata", "example_run", file, package = "magreport")
}

test_that("plot_mag_quality builds a plot", {
  skip_if_not_installed("ggplot2")

  mags <- read_mag_quality(example_path("mag_quality.tsv"))
  p <- plot_mag_quality(mags)

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplot_build(plot_mag_quality(mags, facet_by = "binner")))
  expect_no_error(ggplot2::ggplot_build(plot_mag_quality(mags, label_high = TRUE)))
  expect_error(plot_mag_quality(mags, facet_by = "not_a_column"), "not found")
})

test_that("plot_mag_quality needs the CheckM2 columns", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_mag_quality(data.frame(bin = "a")),
    "completeness"
  )
})

test_that("plot_mapping_rates builds a plot", {
  skip_if_not_installed("ggplot2")

  mapping <- read_mapping_summary(example_path("mapping_summary.tsv"))
  p <- plot_mapping_rates(mapping)

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplot_build(plot_mapping_rates(mapping, warn_below = NULL)))
})

test_that("plot_taxonomy groups unclassified MAGs", {
  skip_if_not_installed("ggplot2")

  mags <- read_mag_quality(example_path("mag_quality.tsv"))

  p <- plot_taxonomy(mags, rank = "phylum")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  # species is empty in the example data
  expect_error(plot_taxonomy(mags, rank = "species"), "No MAG has a taxonomic")

  # genus has both classified MAGs and NA values
  built <- ggplot2::ggplot_build(plot_taxonomy(mags, rank = "genus", fill_by = "quality"))
  expect_gt(nrow(built$data[[1]]), 0)

  expect_error(plot_taxonomy(mags, rank = "kingdom"), "should be one of")
})

test_that("plot_taxonomy explains a missing taxonomy column", {
  skip_if_not_installed("ggplot2")

  mags <- read_mag_quality(example_path("mag_quality.tsv"))
  mags$phylum <- NULL

  expect_error(plot_taxonomy(mags, rank = "phylum"), "gtdbtk_db")
})

test_that("plot_mapping_rates handles samples mapped against several assemblies", {
  skip_if_not_installed("ggplot2")

  mapping <- data.frame(
    assembly = rep(c("T0", "T1"), each = 2),
    sample = rep(c("A", "B"), times = 2),
    read_pairs = 1000,
    overall_alignment_rate_pct = c(90, 20, 15, 85)
  )

  p <- plot_mapping_rates(mapping)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_s3_class(p$facet, "FacetWrap")
})
