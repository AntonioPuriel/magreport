mh_path <- function(file) {
  system.file("extdata", "marine_honour", file, package = "magreport")
}

test_that("read_mimag_quality types the columns and splits the taxonomy", {
  mags <- read_mimag_quality(mh_path("mimag_quality.tsv"))

  expect_equal(nrow(mags), 74)
  expect_equal(mags$bin, mags$genome)
  expect_s3_class(mags$quality, "factor")
  expect_equal(levels(mags$quality), quality_levels())
  expect_type(mags$trna_types, "double")
  expect_true(all(c("phylum", "class", "order", "family", "genus", "species") %in% names(mags)))

  # GTDB-Tk leaves empty ranks as "g__": they become NA
  arch <- mags[mags$genome == "T2.metabat2.344", ]
  expect_equal(arch$domain, "Archaea")
  expect_equal(arch$phylum, "Thermoplasmatota")
  expect_true(is.na(arch$species))

  # Works with the functions written for read_mag_quality()
  counts <- count_quality(mags)
  expect_equal(counts$n[counts$quality == "high"], 6)
})

test_that("split_gtdb_classification handles missing and partial strings", {
  ranks <- magreport:::split_gtdb_classification(c(
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__;f__;g__;s__",
    NA,
    "d__Archaea"
  ))
  expect_equal(ranks$phylum, c("Pseudomonadota", NA, NA))
  expect_equal(ranks$domain, c("Bacteria", NA, "Archaea"))
  expect_true(all(is.na(ranks$order)))
})

test_that("read_hydrocarbon_genes and hydrocarbon_groups", {
  genes <- read_hydrocarbon_genes(mh_path("hydrocarbon_per_mag.tsv"))
  expect_setequal(unique(genes$database), c("CANT-HYD", "HADEG"))
  expect_type(genes$n_genes, "double")

  groups <- hydrocarbon_groups(genes)
  expect_named(groups, c("bin", "group"))
  expect_false(any(duplicated(groups)))
  expect_equal(sum(groups$group == "Any CANT-HYD marker"), 62)

  # AhyA is anaerobic, CYP153 aerobic alkanes, NdoB aromatics
  ahya <- unique(genes$bin[genes$marker == "AhyA"])
  expect_true(all(ahya %in% groups$bin[groups$group == "Anaerobic"]))
  cyp <- unique(genes$bin[genes$marker == "CYP153"])
  expect_true(all(cyp %in% groups$bin[groups$group == "Alkanes (aerobic)"]))
  expect_false("Other" %in% groups$group)

  expect_false("Any CANT-HYD marker" %in% hydrocarbon_groups(genes, any = FALSE)$group)
  expect_error(hydrocarbon_groups(genes[genes$database == "HADEG", ]), "No CANT-HYD")
})

test_that("read_mag_abundance computes relative abundance per sample", {
  ab <- read_mag_abundance(mh_path("mag_abundance.tsv.gz"))
  expect_equal(length(unique(ab$sample)), 48)
  expect_false(any(duplicated(ab[c("bin", "sample")])))
  sums <- tapply(ab$rel_abundance, ab$sample, sum)
  expect_equal(as.vector(sums), rep(100, 48))

  two <- unique(ab$bin)[1:2]
  sub <- read_mag_abundance(mh_path("mag_abundance.tsv.gz"), bins = two)
  expect_setequal(unique(sub$bin), two)
  expect_equal(as.vector(tapply(sub$rel_abundance, sub$sample, sum)), rep(100, 48))

  expect_warning(
    read_mag_abundance(mh_path("mag_abundance.tsv.gz"), bins = c(two, "not_a_bin")),
    "not found"
  )
})

test_that("parse_sample_names splits the design", {
  s <- parse_sample_names(c("T0-BB-C-1", "T7-BB-O-3", "T0-BB-C-1"),
                          into = c("time", "site", "condition", "replicate"))
  expect_equal(nrow(s), 2)
  expect_equal(s$condition, c("C", "O"))
  expect_error(parse_sample_names("T0-BB", into = c("a", "b", "c")), "3 parts")
})

test_that("plot_abundance_timeline summarises groups over time", {
  skip_if_not_installed("ggplot2")

  ab <- read_mag_abundance(mh_path("mag_abundance.tsv.gz"))
  genes <- read_hydrocarbon_genes(mh_path("hydrocarbon_per_mag.tsv"))
  samples <- parse_sample_names(ab$sample, into = c("time", "site", "condition", "replicate"))
  samples$day <- c(T0 = 2, T1 = 17, T2 = 33, T3 = 48, T4 = 63, T5 = 93, T6 = 124, T7 = 185)[samples$time]

  p <- plot_abundance_timeline(ab, samples, groups = hydrocarbon_groups(genes),
                               time = "day", condition = "condition")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  # 5 groups x 8 times x 2 conditions, 3 replicates each
  expect_equal(nrow(p$data), 5 * 8 * 2)
  expect_true(all(p$data$n == 3))

  # The sum over all MAGs is 100 % in every sample
  total <- plot_abundance_timeline(ab, samples, time = "time")
  expect_equal(total$data$mean, rep(100, 8))
  expect_equal(total$data$sd, rep(0, 8), tolerance = 1e-8)
  expect_no_error(ggplot2::ggplot_build(total))

  labelled <- plot_abundance_timeline(ab, samples, time = "day", condition = "condition",
                                      condition_labels = c(C = "Low impact", O = "High impact"),
                                      condition_title = "Site")
  built <- ggplot2::ggplot_build(labelled)
  colour <- built$plot$scales$get_scales("colour")
  expect_equal(as.character(colour$get_labels()), c("Low impact", "High impact"))
  expect_equal(colour$name, "Site")
  expect_error(
    plot_abundance_timeline(ab, samples, condition = "condition", condition_labels = c(C = "x")),
    "no label for: O"
  )

  expect_error(plot_abundance_timeline(ab, samples[-1, ]), "missing from samples")
  expect_error(plot_abundance_timeline(ab, samples, time = "nope"), "not found")
})

test_that("plot_hydrocarbon_genes builds a heatmap", {
  skip_if_not_installed("ggplot2")

  genes <- read_hydrocarbon_genes(mh_path("hydrocarbon_per_mag.tsv"))
  mags <- read_mimag_quality(mh_path("mimag_quality.tsv"))

  p <- plot_hydrocarbon_genes(genes, mags = mags)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(length(unique(p$data$bin)), 40)

  hadeg <- plot_hydrocarbon_genes(genes, database = "HADEG", top = 10)
  expect_equal(length(unique(hadeg$data$bin)), 10)
  expect_no_error(ggplot2::ggplot_build(hadeg))

  expect_error(plot_hydrocarbon_genes(genes, mags = mags, rank = "nope"), "not found")
})
