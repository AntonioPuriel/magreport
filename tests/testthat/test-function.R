example_path <- function(file) {
  system.file("extdata", "example_run", file, package = "magreport")
}

test_that("read_ko_table and read_annotation_stats return typed columns", {
  kos <- read_ko_table(example_path("ko_per_mag.tsv"))
  stats <- read_annotation_stats(example_path("annotation_stats.tsv"))

  expect_true(all(c("bin", "ko", "n_genes") %in% names(kos)))
  expect_type(kos$n_genes, "double")

  expect_true("annotated_pct" %in% names(stats))
  expect_true(all(stats$annotated_pct > 0 & stats$annotated_pct <= 100))
  expect_equal(stats$annotated_pct, 100 * stats$n_annotated / stats$n_genes)
})

test_that("hydrocarbon_kos is a usable marker set", {
  markers <- hydrocarbon_kos()

  expect_named(markers, c("pathway", "ko"))
  expect_gt(nrow(markers), 30)
  expect_false(any(duplicated(markers$ko)))
  expect_true(all(grepl("^K[0-9]{5}$", markers$ko)))
})

test_that("plot_pathway_presence builds a plot", {
  skip_if_not_installed("ggplot2")

  kos <- read_ko_table(example_path("ko_per_mag.tsv"))
  mags <- read_mag_quality(example_path("mag_quality.tsv"))

  expect_s3_class(plot_pathway_presence(kos), "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot_pathway_presence(kos, mags = mags)))
  expect_error(plot_pathway_presence(kos, min_kos = 1e6), "min_kos")

  other <- data.frame(pathway = "Nothing", ko = "K99999")
  expect_error(plot_pathway_presence(kos, markers = other), "None of the marker genes")
  expect_error(plot_pathway_presence(kos, markers = data.frame(a = 1)), "pathway and ko")
})

test_that("read_pipeline picks up the functional tables", {
  run <- read_pipeline(system.file("extdata", "example_run", package = "magreport"))

  expect_s3_class(run$kos, "data.frame")
  expect_s3_class(run$annotation, "data.frame")
})
