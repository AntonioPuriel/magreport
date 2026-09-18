example_path <- function(file) {
  system.file("extdata", "example_run", file, package = "magreport")
}

test_that("read_bin_summary reads every bin", {
  bins <- read_bin_summary(example_path("bin_summary.tsv"))

  expect_true(all(c("binner", "bin", "n_contigs", "total_length") %in% names(bins)))
  expect_type(bins$total_length, "double")
  expect_true(all(bins$n_contigs > 0))
})

test_that("read_mapping_summary reads one row per sample", {
  mapping <- read_mapping_summary(example_path("mapping_summary.tsv"))

  expect_equal(nrow(mapping), 2)
  expect_type(mapping$read_pairs, "double")
  expect_true(all(mapping$overall_alignment_rate_pct >= 0))
  expect_true(all(mapping$overall_alignment_rate_pct <= 100))
})

test_that("read_depth reshapes the depth table", {
  wide <- read_depth(example_path("coassembly.depth.txt"), long = FALSE)
  long <- read_depth(example_path("coassembly.depth.txt"))

  expect_true("totalAvgDepth" %in% names(wide))
  expect_identical(names(long), c("contig", "contig_length", "sample", "depth", "depth_var"))
  expect_setequal(unique(long$sample), c("minigut1", "minigut2"))
  expect_equal(nrow(long), nrow(wide) * 2)
  expect_false(anyNA(long$depth))
})
