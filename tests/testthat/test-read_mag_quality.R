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

test_that("N/A from GTDB-Tk is read as a missing value", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(
    c(
      paste("assembly", "binner", "bin", "n_contigs", "total_length", "n50",
            "gc_pct", "completeness", "contamination", "quality", "domain",
            "phylum", "closest_ani", sep = "\t"),
      paste("asm", "dastool", "bin.1", "10", "1000", "500", "45", "95", "1",
            "high", "Bacteria", "Pseudomonadota", "98.02", sep = "\t"),
      paste("asm", "dastool", "bin.2", "5", "500", "200", "40", "30", "0",
            "low", "NA", "NA", "N/A", sep = "\t")
    ),
    path
  )

  mags <- read_mag_quality(path)

  expect_true(is.na(mags$closest_ani[2]))
  expect_true(is.na(mags$phylum[2]))
  expect_equal(mags$closest_ani[1], 98.02)
  expect_type(mags$closest_ani, "double")
})
