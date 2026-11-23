extdata <- function(...) system.file("extdata", ..., package = "magreport")

# Two copies of the Marine Honour GUNC and MIMAG tables: the same bin names
# in two studies, as can happen in a catalog
local_two_studies <- function(env = parent.frame()) {
  dirs <- vapply(c("a", "b"), function(study) {
    dir <- withr::local_tempdir(.local_envir = env)
    dir.create(file.path(dir, "05_quality", "gunc"), recursive = TRUE)
    file.copy(list.files(extdata("marine_honour", "gunc"), full.names = TRUE),
              file.path(dir, "05_quality", "gunc"))
    file.copy(extdata("marine_honour", "mimag_quality.tsv"), file.path(dir, "05_quality"))
    dir
  }, character(1))
  # Study b: a different CheckM2 contamination for every MAG
  path <- file.path(dirs[["b"]], "05_quality", "mimag_quality.tsv")
  x <- utils::read.delim(path, stringsAsFactors = FALSE)
  x$contamination <- x$contamination + 1
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  dirs
}

test_that("read_pipeline stacks several runs with a study column", {
  runs <- read_pipeline(c(gut = extdata("example_run"), sediment = extdata("marine_honour")))
  gut <- read_pipeline(extdata("example_run"))
  sediment <- read_pipeline(extdata("marine_honour"))

  expect_identical(attr(runs, "studies"), c("gut", "sediment"))
  expect_identical(names(runs), names(gut))
  expect_identical(names(runs$mags)[1], "study")
  expect_equal(nrow(runs$mags), nrow(gut$mags))
  expect_equal(nrow(runs$mimag), nrow(sediment$mimag))
  expect_identical(unique(runs$mimag$study), "sediment")
  expect_equal(runs$mags[, -1], gut$mags, ignore_attr = TRUE)
  expect_s3_class(runs$mimag$quality, "factor")
  expect_identical(names(attr(runs, "files")), c("gut", "sediment"))
  expect_null(runs$singlem)
})

test_that("read_pipeline fills the columns that only some studies have", {
  dirs <- local_two_studies()
  path <- file.path(dirs[["a"]], "05_quality", "mimag_quality.tsv")
  x <- utils::read.delim(path, stringsAsFactors = FALSE)
  x$extra <- "only in a"
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  runs <- read_pipeline(dirs)
  expect_identical(attr(runs, "studies"), c("a", "b"))
  expect_true(all(runs$mimag$extra[runs$mimag$study == "a"] == "only in a"))
  expect_true(all(is.na(runs$mimag$extra[runs$mimag$study == "b"])))
})

test_that("read_pipeline needs unique study names", {
  dir <- extdata("example_run")
  expect_error(read_pipeline(c(dir, dir)), "unique")
  expect_error(read_pipeline(c(a = dir, a = dir)), "unique")
  expect_error(read_pipeline(stats::setNames(c(dir, dir), c("a", ""))), "study name")
  expect_error(read_pipeline(character()), "No results directory")
  # A single named directory also gets the study column
  one <- read_pipeline(c(gut = dir))
  expect_identical(unique(one$mags$study), "gut")
})

test_that("joins by bin stay within each study", {
  runs <- read_pipeline(local_two_studies())
  expect_equal(nrow(runs$gunc), 2 * 1161)

  p <- plot_gunc(runs$gunc, mags = runs$mimag)
  data <- p$data
  expect_equal(nrow(data), 2 * 281)
  # Same bin, each study with its own contamination
  bin <- runs$mimag$bin[1]
  expect_equal(
    data$x[data$study == "b" & data$bin == bin] - data$x[data$study == "a" & data$bin == bin],
    1
  )
})

test_that("mag_recovery uses the MAGs of each study", {
  profile <- data.frame(
    study = rep(c("a", "b"), each = 2),
    rank = "genus",
    taxon = c("g__Alcanivorax", "g__Marinobacter", "g__Alcanivorax", "g__Marinobacter"),
    sample = rep(c("s1", "s2"), each = 2),
    rel_abundance = c(60, 40, 30, 70),
    stringsAsFactors = FALSE
  )
  mags <- data.frame(
    study = c("a", "b"),
    bin = "m1",
    quality = "high",
    genus = c("Alcanivorax", "Marinobacter"),
    stringsAsFactors = FALSE
  )
  r <- mag_recovery(profile, mags, ranks = "genus")
  expect_identical(r$study, c("a", "b"))
  expect_equal(r$recovered, c(60, 70))

  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_mag_recovery(r), "ggplot")
})

test_that("summarise_studies gives one row per study", {
  runs <- read_pipeline(c(gut = extdata("example_run"), sediment = extdata("marine_honour")))
  s <- summarise_studies(runs)

  expect_identical(s$study, c("gut", "sediment"))
  expect_equal(s$samples[1], 2)
  expect_equal(sum(unlist(s[2, c("high", "medium", "low")])), 281)
  expect_true(is.na(s$gunc_fail[1]))
  expect_true(is.na(s$bins[2]))
  expect_error(summarise_studies(read_pipeline(extdata("example_run"))), "single run")
})

test_that("mag_report builds a report of several studies", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  out <- withr::local_tempfile(fileext = ".html")
  dirs <- c(gut = extdata("example_run"), sediment = extdata("marine_honour"))

  mag_report(dirs, out)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "Studies. MAGs per quality")
  expect_match(html, "sediment")
})
