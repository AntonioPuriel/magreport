contigs_path <- function(...) {
  system.file("extdata", "example_contigs", ..., package = "magreport")
}

example_samples <- function(profile) {
  parse_sample_names(profile$sample, into = c("time", "site", "condition", "replicate"))
}

test_that("read_taxonomy_profile adds up to 100 % per sample and rank", {
  cp <- read_taxonomy_profile(contigs_path("profiles", "taxonomy_per_sample.tsv.gz"))

  expect_true(all(c("rank", "taxon", "sample", "n_genes", "tpm", "rel_abundance") %in% names(cp)))
  expect_identical(attr(cp, "source"), "CAT")
  expect_setequal(unique(cp$rank), taxonomic_ranks())
  totals <- tapply(cp$rel_abundance, paste(cp$sample, cp$rank), sum)
  expect_true(all(abs(totals - 100) < 1e-8))
  expect_true("Unclassified" %in% cp$taxon)
})

test_that("read_taxonomy_profile rejects unknown ranks", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("rank\ttaxon\tsample\tn_genes\ttpm", "kingdom\tk__X\ts1\t1\t10"), path)
  expect_error(read_taxonomy_profile(path), "Unexpected rank")
})

test_that("read_singlem builds every rank from the profiles of each sample", {
  files <- list.files(contigs_path("singlem"), full.names = TRUE)
  sm <- read_singlem(files)

  expect_identical(attr(sm, "source"), "SingleM")
  expect_identical(names(sm), c("rank", "taxon", "sample", "coverage", "rel_abundance"))
  expect_length(unique(sm$sample), length(files))
  totals <- tapply(sm$rel_abundance, paste(sm$sample, sm$rank), sum)
  expect_true(all(abs(totals - 100) < 1e-8))

  # The coverage of a taxon includes that of its descendants
  one <- sm[sm$sample == "T0-BB-C-1", ]
  raw <- read.delim(contigs_path("singlem", "T0-BB-C-1.singlem.tsv"))
  expect_equal(
    one$coverage[one$rank == "domain" & one$taxon == "d__Archaea"],
    sum(raw$coverage[grepl("d__Archaea", raw$taxonomy)])
  )
  # Reads classified only to Root are unclassified at every rank
  expect_gt(one$coverage[one$rank == "domain" & one$taxon == "Unclassified"], 0)
})

test_that("read_singlem reads a taxon-by-sample table", {
  path <- file.path(withr::local_tempdir(), "singlem_phylum.tsv")
  writeLines(c(
    "taxonomy\tT0-BB-C-1\tT0-BB-O-1",
    "Root; d__Bacteria; p__Pseudomonadota\t50\t30",
    "Root; d__Bacteria; p__Bacteroidota\t20\t40",
    "Root; d__Bacteria\t10\t10"
  ), path)
  sm <- read_singlem(path)

  expect_setequal(unique(sm$sample), c("T0-BB-C-1", "T0-BB-O-1"))
  expect_identical(unique(sm$rank), "phylum")
  totals <- tapply(sm$rel_abundance, sm$sample, sum)
  expect_true(all(abs(totals - 100) < 1e-8))
  # 10 % stops at the domain and 20 % is missing: 30 % unclassified
  un <- sm[sm$taxon == "Unclassified" & sm$sample == "T0-BB-C-1", "rel_abundance"]
  expect_equal(un, 30)
})

test_that("read_feature_profile and read_profile_stats read the functional profiles", {
  stats <- read_profile_stats(contigs_path("profiles", "profile_stats.tsv"))
  expect_true(all(grepl("^fraction_", setdiff(names(stats), c("sample", "total_tpm")))))
  expect_true(all(stats$fraction_cat_domain > 0 & stats$fraction_cat_domain <= 1))

  markers <- read_feature_profile(contigs_path("profiles", "canthyd_per_sample.tsv.gz"),
                                  total_tpm = stats)
  expect_identical(markers$feature, markers$marker)
  expect_equal(markers$rel_abundance, 100 * markers$tpm / 1e6)

  kos <- read_feature_profile(contigs_path("profiles", "ko_per_sample.tsv.gz"))
  expect_identical(kos$feature, kos$ko)
  expect_error(
    read_feature_profile(contigs_path("profiles", "ko_per_sample.tsv.gz"),
                         total_tpm = c(other = 1)),
    "no value for sample"
  )
})

test_that("read_catalog_stats reads one numeric row", {
  catalog <- read_catalog_stats(contigs_path("catalog", "gene_catalog.stats.tsv"))
  expect_equal(nrow(catalog), 1)
  expect_type(catalog$catalog_genes, "double")
})

test_that("read_pipeline finds a run without MAGs", {
  run <- read_pipeline(contigs_path())

  expect_null(run$mags)
  for (name in c("catalog", "profile_stats", "contig_taxonomy", "contig_kos",
                 "contig_canthyd", "singlem")) {
    expect_false(is.null(run[[name]]), info = name)
  }
  expect_null(run$contig_hadeg)
  expect_length(attr(run, "files")$singlem, 12)
})

test_that("plot_community draws CAT and SingleM profiles", {
  skip_if_not_installed("ggplot2")
  cp <- read_taxonomy_profile(contigs_path("profiles", "taxonomy_per_sample.tsv.gz"))
  sm <- read_singlem(list.files(contigs_path("singlem"), full.names = TRUE))
  samples <- example_samples(cp)

  p <- plot_community(cp, rank = "genus", top = 3)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_setequal(levels(p$data$taxon),
                  c("g__Alcanivorax", "g__Polaribacter", "g__Marinobacter", "Other", "Unclassified"))

  expect_no_error(ggplot2::ggplot_build(
    plot_community(sm, rank = "phylum", samples = samples, facet_by = "condition")
  ))
  expect_error(plot_community(cp, facet_by = "condition"), "needs the samples")
  expect_error(plot_community(cp, rank = "kingdom"), "should be one of")
  expect_error(plot_community(data.frame(a = 1)), "must come from")
})

test_that("plot_profile_comparison compares the classified parts", {
  skip_if_not_installed("ggplot2")
  cp <- read_taxonomy_profile(contigs_path("profiles", "taxonomy_per_sample.tsv.gz"))
  sm <- read_singlem(list.files(contigs_path("singlem"), full.names = TRUE))

  p <- plot_profile_comparison(cp, sm, rank = "phylum")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  # Each profile is rescaled to its classified part: the means add up to 100
  expect_equal(sum(p$data$mean_x), 100)
  expect_equal(sum(p$data$mean_y), 100)
  expect_match(p$labels$subtitle, "Spearman")

  other <- sm
  other$sample <- paste0("x", other$sample)
  expect_error(plot_profile_comparison(cp, other), "no sample in common")
})

test_that("plot_annotation_fraction and plot_feature_heatmap build plots", {
  skip_if_not_installed("ggplot2")
  stats <- read_profile_stats(contigs_path("profiles", "profile_stats.tsv"))
  expect_no_error(ggplot2::ggplot_build(plot_annotation_fraction(stats)))
  expect_error(plot_annotation_fraction(data.frame(sample = "a")), "must come from")

  kos <- read_feature_profile(contigs_path("profiles", "ko_per_sample.tsv.gz"))
  samples <- example_samples(kos)
  p <- plot_feature_heatmap(kos, features = hydrocarbon_kos(), samples = samples,
                            facet_by = "condition")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_true("pathway" %in% names(p$data))
  expect_no_error(ggplot2::ggplot_build(plot_feature_heatmap(kos, top = 5)))
  expect_error(plot_feature_heatmap(kos, features = "K99999"), "None of the features")
})

test_that("plot_abundance_timeline follows features and taxa", {
  skip_if_not_installed("ggplot2")
  markers <- read_feature_profile(contigs_path("profiles", "canthyd_per_sample.tsv.gz"))
  samples <- example_samples(markers)

  p <- plot_abundance_timeline(markers, samples, time = "time", condition = "condition")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_setequal(as.character(unique(p$data$group)), unique(markers$feature))

  kos <- read_feature_profile(contigs_path("profiles", "ko_per_sample.tsv.gz"))
  p <- plot_abundance_timeline(kos, samples, groups = hydrocarbon_kos(),
                               time = "time", condition = "condition")
  expect_true("Alkanes (aerobic)" %in% p$data$group)

  cp <- read_taxonomy_profile(contigs_path("profiles", "taxonomy_per_sample.tsv.gz"))
  expect_error(plot_abundance_timeline(cp, samples, time = "time"), "several ranks")
  genus <- cp[cp$rank == "genus" & cp$taxon != "Unclassified", ]
  expect_no_error(ggplot2::ggplot_build(
    plot_abundance_timeline(genus, samples, time = "time", condition = "condition")
  ))
})
