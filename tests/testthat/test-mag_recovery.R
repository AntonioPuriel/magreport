toy_profile <- function() {
  data.frame(
    rank = rep(c("phylum", "genus"), each = 6),
    taxon = c("p__Bacteroidota", "p__Pseudomonadota", "Unclassified",
              "p__Bacteroidota", "p__Pseudomonadota", "Unclassified",
              "g__Polaribacter", "g__Alcanivorax", "Unclassified",
              "g__Polaribacter", "g__Marinobacter", "Unclassified"),
    sample = rep(rep(c("A", "B"), each = 3), 2),
    rel_abundance = c(50, 40, 10,  20, 70, 10,  30, 20, 50,  10, 60, 30),
    stringsAsFactors = FALSE
  )
}

toy_mags <- function() {
  data.frame(
    bin = c("m1", "m2", "m3"),
    quality = factor(c("high", "medium", "low"), levels = quality_levels()),
    phylum = c("Bacteroidota", "Pseudomonadota", "Pseudomonadota"),
    genus = c("Polaribacter", NA, "Alcanivorax"),
    stringsAsFactors = FALSE
  )
}

test_that("mag_recovery sums the share of taxa with a MAG", {
  r <- mag_recovery(toy_profile(), toy_mags(), ranks = c("phylum", "genus"))

  expect_identical(levels(r$rank), c("phylum", "genus"))
  get <- function(s, k, col) r[[col]][r$sample == s & r$rank == k]
  expect_equal(get("A", "phylum", "recovered"), 90)
  expect_equal(get("A", "phylum", "classified"), 90)
  # The Alcanivorax MAG is of low quality: left out by default
  expect_equal(get("A", "genus", "recovered"), 30)
  expect_equal(get("A", "genus", "classified"), 50)
  expect_equal(get("A", "genus", "recovered_of_classified"), 60)
  expect_equal(get("B", "genus", "recovered"), 10)
  expect_equal(get("B", "genus", "n_taxa_with_mag"), 1L)

  every <- mag_recovery(toy_profile(), toy_mags(), ranks = c("phylum", "genus"),
                        min_quality = "low")
  expect_equal(every$recovered[every$sample == "A" & every$rank == "genus"], 50)
})

test_that("mag_recovery checks its inputs", {
  expect_error(mag_recovery(data.frame(x = 1), toy_mags()), "read_singlem")
  expect_error(mag_recovery(toy_profile(), data.frame(bin = "a")), "no taxonomic rank")
  expect_error(mag_recovery(toy_profile(), toy_mags(), ranks = "order"), "no taxonomic rank")
  mags <- toy_mags()
  mags$order <- "x"
  expect_error(mag_recovery(toy_profile(), mags, ranks = "order"), "none of the ranks")
})

test_that("mag_recovery works on the example SingleM profile", {
  singlem <- read_singlem(list.files(
    system.file("extdata", "example_contigs", "singlem", package = "magreport"),
    full.names = TRUE
  ))
  mimag <- read_mimag_quality(system.file("extdata", "marine_honour", "mimag_quality.tsv",
                                          package = "magreport"))
  r <- mag_recovery(singlem, mimag)

  expect_equal(nrow(r), length(unique(singlem$sample)) * 5)
  expect_true(all(r$recovered <= r$classified + 1e-9))
  expect_true(all(r$classified <= 100 + 1e-9))
  expect_identical(attr(r, "source"), "SingleM")

  skip_if_not_installed("ggplot2")
  p <- plot_mag_recovery(r)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "SingleM")
  expect_error(plot_mag_recovery(data.frame(x = 1)), "mag_recovery")
})
