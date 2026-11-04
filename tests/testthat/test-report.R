example_dir <- function() {
  system.file("extdata", "example_run", package = "magreport")
}

test_that("read_pipeline finds every table", {
  run <- read_pipeline(example_dir())

  expect_named(run, c(
    "mags", "bins", "mapping", "depth", "assembly", "kos", "annotation",
    "mimag", "hydrocarbon", "abundance"
  ))
  expect_s3_class(run$mags, "data.frame")
  expect_s3_class(run$bins, "data.frame")
  expect_equal(nrow(run$mapping), 2)
  files <- attr(run, "files")
  older <- c("mags", "bins", "mapping", "depth", "assembly", "kos", "annotation")
  expect_true(all(vapply(files[older], Negate(is.null), logical(1))))

  # Tables of mag-pipeline >= 1.8.0, not in this older example
  expect_null(run$mimag)
  expect_null(run$abundance)
})

test_that("read_pipeline returns NULL for missing tables", {
  tmp <- withr::local_tempdir()
  file.copy(file.path(example_dir(), "mag_quality.tsv"), tmp)

  run <- read_pipeline(tmp)

  expect_s3_class(run$mags, "data.frame")
  expect_null(run$mapping)
  expect_null(run$depth)

  expect_error(read_pipeline(file.path(tmp, "nope")), "Directory not found")
})

test_that("read_pipeline reads the MIMAG, hydrocarbon and abundance tables", {
  run <- read_pipeline(system.file("extdata", "marine_honour", package = "magreport"))

  expect_s3_class(run$mimag, "data.frame")
  expect_s3_class(run$hydrocarbon, "data.frame")
  # abundance is restricted to the genomes of the MIMAG table
  expect_setequal(unique(run$abundance$bin), run$mimag$bin)
  sums <- tapply(run$abundance$rel_abundance, run$abundance$sample, sum)
  expect_equal(as.vector(sums), rep(100, length(sums)))
})

test_that("mag_report writes an HTML file", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("ggplot2")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")

  out <- file.path(withr::local_tempdir(), "report.html")
  mag_report(example_dir(), out)

  expect_true(file.exists(out))
  expect_gt(file.size(out), 10000)

  html <- readLines(out, warn = FALSE)
  expect_true(any(grepl("MAG quality", html)))
})

test_that("mag_report fails clearly on an empty directory", {
  skip_if_not_installed("rmarkdown")

  empty <- withr::local_tempdir()
  expect_error(
    mag_report(empty, file.path(empty, "report.html")),
    "No mag-pipeline tables found"
  )
})
