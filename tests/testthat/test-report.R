example_dir <- function() {
  system.file("extdata", "example_run", package = "magreport")
}

test_that("read_pipeline finds every table", {
  run <- read_pipeline(example_dir())

  expect_named(run, c("mags", "bins", "mapping", "depth", "assembly"))
  expect_s3_class(run$mags, "data.frame")
  expect_s3_class(run$bins, "data.frame")
  expect_equal(nrow(run$mapping), 2)
  expect_true(all(vapply(attr(run, "files"), Negate(is.null), logical(1))))
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
