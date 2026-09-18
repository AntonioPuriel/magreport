example_path <- function(file) {
  system.file("extdata", "example_run", file, package = "magreport")
}

test_that("read_mag_quality returns typed columns", {
  mags <- read_mag_quality(example_path("mag_quality.tsv"))

  expect_s3_class(mags, "data.frame")
  expect_true(all(c("assembly", "binner", "bin", "completeness") %in% names(mags)))
  expect_type(mags$completeness, "double")
  expect_s3_class(mags$quality, "factor")
  expect_identical(levels(mags$quality), quality_levels())
  expect_gt(nrow(mags), 0)
})

test_that("missing files and malformed tables are reported", {
  expect_error(read_mag_quality("does_not_exist.tsv"), "File not found")

  bad <- tempfile(fileext = ".tsv")
  on.exit(unlink(bad), add = TRUE)
  write.table(
    data.frame(assembly = "a", binner = "b"),
    bad,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  expect_error(read_mag_quality(bad), "missing column")
})

test_that("count_quality keeps empty categories", {
  mags <- read_mag_quality(example_path("mag_quality.tsv"))
  counts <- count_quality(mags)

  expect_identical(as.character(counts$quality), quality_levels())
  expect_equal(sum(counts$n), nrow(mags))

  by_binner <- count_quality(mags, by = "binner")
  expect_true("binner" %in% names(by_binner))
  expect_equal(sum(by_binner$n), nrow(mags))

  expect_error(count_quality(mags, by = "not_a_column"), "not found")
})
