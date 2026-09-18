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

test_that("read_depth copes with different depth headers", {
  write_depth <- function(header) {
    path <- withr::local_tempfile(fileext = ".depth.txt", .local_envir = parent.frame())
    writeLines(c(header, "k141_1\t1325\t3.0\t1.0\t0.5\t5.0\t2.5"), path)
    path
  }

  headers <- c(
    "contigName\tcontigLen\ttotalAvgDepth\tA.sorted.bam\tA.sorted.bam-var\tB.sorted.bam\tB.sorted.bam-var",
    "contigName\tcontigLen\ttotalAvgDepth\tA\tA-var\tB\tB-var",
    "contigName\tcontigLen\ttotalAvgDepth\tA.bam\tA.bam.var\tB.bam\tB.bam.var",
    "contigName\tcontigLen\ttotalAvgDepth\tA.sorted\tA.sorted-var\tB.sorted\tB.sorted-var"
  )

  for (header in headers) {
    depth <- read_depth(write_depth(header))
    expect_setequal(unique(depth$sample), c("A", "B"))
    expect_equal(depth$depth, c(1, 5))
    expect_equal(depth$depth_var, c(0.5, 2.5))
  }
})

test_that("read_depth reports a file without sample columns", {
  path <- withr::local_tempfile(fileext = ".depth.txt")
  writeLines(c("contigName\tcontigLen\ttotalAvgDepth", "k141_1\t1325\t3.0"), path)

  expect_error(read_depth(path), "No per-sample depth columns")
})
