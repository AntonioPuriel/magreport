extdata <- function(...) system.file("extdata", ..., package = "magreport")

read_out <- function(dir, file) {
  utils::read.delim(file.path(dir, file), check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("export_tables writes the community profile as matrices", {
  dir <- withr::local_tempdir()
  files <- export_tables(extdata("example_contigs"), dir, excel = FALSE, phyloseq = FALSE)

  expect_s3_class(files, "data.frame")
  expect_named(files, c("file", "rows", "columns", "description"))
  expect_true(all(file.exists(file.path(dir, files$file))))
  expect_true("README.txt" %in% files$file)

  samples <- read_out(dir, "samples.tsv")
  expect_equal(nrow(samples), 12)
  expect_true(all(c("sample", "total_tpm", "singlem_coverage") %in% names(samples)))

  # Both measures of every rank; relative columns add up to 100
  for (source in c("cat", "singlem")) {
    raw <- if (source == "cat") "tpm" else "coverage"
    for (rank in c("domain", "phylum", "genus", "species")) {
      expect_true(file.exists(file.path(dir, "matrices", sprintf("%s_%s_%s.tsv", source, rank, raw))))
      rel <- read_out(dir, sprintf("matrices/%s_%s_relative.tsv", source, rank))
      expect_identical(names(rel)[1], "taxon")
      expect_equal(unname(colSums(rel[-1])), rep(100, 12), tolerance = 1e-6)
    }
  }

  # The matrix holds the same values as the long table
  run <- read_pipeline(extdata("example_contigs"))
  long <- run$singlem[run$singlem$rank == "phylum", ]
  wide <- read_out(dir, "matrices/singlem_phylum_coverage.tsv")
  one <- long[long$sample == "T1-BB-O-1" & long$taxon == "p__Pseudomonadota", "coverage"]
  expect_equal(wide[wide$taxon == "p__Pseudomonadota", "T1-BB-O-1"], one)

  expect_true(file.exists(file.path(dir, "matrices", "ko_tpm.tsv")))
  expect_true(file.exists(file.path(dir, "matrices", "canthyd_relative.tsv")))
  expect_true(file.exists(file.path(dir, "long", "contig_taxonomy.tsv")))
  expect_false(file.exists(file.path(dir, "mags.tsv")))

  readme <- readLines(file.path(dir, "README.txt"))
  expect_true(any(grepl("singlem_phylum_relative.tsv", readme, fixed = TRUE)))
  expect_true(any(grepl("no phyloseq object", readme)))
})

test_that("export_tables summarises the MAGs and their abundance", {
  dir <- withr::local_tempdir()
  export_tables(extdata("marine_honour"), dir, excel = FALSE, phyloseq = FALSE)

  mags <- read_out(dir, "mags.tsv")
  expect_equal(nrow(mags), 281)
  expect_identical(names(mags)[1], "bin")
  expect_false("genome" %in% names(mags))
  expect_true(all(c("quality", "phylum", "abundance_mean_pct", "n_samples_with_reads",
                    "n_canthyd_markers", "canthyd_markers") %in% names(mags)))
  expect_true(all(mags$n_canthyd_markers >= 0))

  depth <- read_out(dir, "matrices/mag_depth.tsv")
  relative <- read_out(dir, "matrices/mag_relative.tsv")
  expect_equal(dim(depth), c(281, 49))
  expect_setequal(depth$bin, mags$bin)
  sums <- colSums(relative[-1])
  expect_equal(unname(sums[sums > 0]), rep(100, sum(sums > 0)), tolerance = 1e-6)

  hc <- read_out(dir, "matrices/mag_hydrocarbon_genes.tsv")
  expect_true(all(grepl("^(CANT-HYD|HADEG):", hc$marker)))

  # Coassembly by group: MAGs have reads in several samples, no warning
  expect_false(any(grepl("one sample only", readLines(file.path(dir, "README.txt")))))
})

test_that("export_tables warns when every MAG has reads in one sample only", {
  run <- read_pipeline(extdata("marine_honour"))
  ab <- run$abundance
  first <- !duplicated(ab$bin)
  ab$mean_depth[!first] <- 0
  run$abundance <- ab

  dir <- withr::local_tempdir()
  export_tables(run, dir, excel = FALSE, phyloseq = FALSE)
  readme <- readLines(file.path(dir, "README.txt"))
  expect_true(any(grepl("one sample only", readme)))
})

test_that("export_tables names samples by study with several runs", {
  run <- read_pipeline(c(mh = extdata("marine_honour"), ej = extdata("example_contigs")),
                       depth = FALSE)
  dir <- withr::local_tempdir()
  export_tables(run, dir, excel = FALSE, phyloseq = FALSE)

  samples <- read_out(dir, "samples.tsv")
  expect_identical(names(samples)[1:2], c("study", "sample"))
  expect_setequal(unique(samples$study), c("mh", "ej"))

  sm <- read_out(dir, "matrices/singlem_phylum_relative.tsv")
  expect_true(all(startsWith(names(sm)[-1], "ej__")))
  depth <- read_out(dir, "matrices/mag_depth.tsv")
  expect_true(all(startsWith(names(depth)[-1], "mh__")))
  expect_true(all(startsWith(depth$bin, "mh__")))
})

test_that("export_tables writes the methods of a run with a manifest", {
  dir <- withr::local_tempdir()
  file.copy(extdata("example_run", "mag_quality.tsv"), dir)
  file.copy(extdata("example_run", "mapping_summary.tsv"), dir)
  manifest <- list(
    schema_version = "1.0",
    pipeline = list(name = "AntonioPuriel/mag-pipeline", version = "2.1.0"),
    run = list(name = "test_run", success = TRUE),
    databases = list(gtdbtk = list(path = "/db/gtdbtk", modified = "2026-05-02")),
    params = list(workflow = "mags", min_contig_length = 500),
    tables = list(mag_quality = I("mag_quality.tsv"), mapping_summary = I("mapping_summary.tsv"))
  )
  jsonlite::write_json(manifest, file.path(dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)

  out <- withr::local_tempdir()
  export_tables(dir, out, excel = FALSE, phyloseq = FALSE)

  params <- read_out(out, "methods/params.tsv")
  expect_true("min_contig_length" %in% params$parameter)
  databases <- read_out(out, "methods/databases.tsv")
  expect_identical(databases$database, "gtdbtk")
  expect_true(any(grepl("mag-pipeline 2.1.0", readLines(file.path(out, "README.txt")))))
})

test_that("export_tables writes an Excel file with one sheet per table", {
  skip_if_not_installed("writexl")
  skip_if_not_installed("readxl")

  dir <- withr::local_tempdir()
  files <- export_tables(extdata("example_contigs"), dir, phyloseq = FALSE)

  xlsx <- file.path(dir, "magreport_tables.xlsx")
  expect_true(file.exists(xlsx))
  sheets <- readxl::excel_sheets(xlsx)
  expect_identical(sheets[1], "README")
  expect_true(all(c("samples", "singlem_phylum_relative", "cat_species_tpm") %in% sheets))
  expect_true(all(nchar(sheets) <= 31))
  expect_false(any(startsWith(sheets, "contig_")))

  phylum <- readxl::read_excel(xlsx, sheet = "singlem_phylum_relative")
  expect_equal(nrow(phylum), nrow(read_out(dir, "matrices/singlem_phylum_relative.tsv")))
})

test_that("export_tables builds phyloseq objects", {
  skip_if_not_installed("phyloseq")

  dir <- withr::local_tempdir()
  export_tables(extdata("example_contigs"), dir, excel = FALSE)
  sm <- readRDS(file.path(dir, "phyloseq", "singlem.rds"))
  expect_s4_class(sm, "phyloseq")
  expect_equal(phyloseq::nsamples(sm), 12)
  expect_true("genus" %in% phyloseq::rank_names(sm))
  expect_true("total_tpm" %in% phyloseq::sample_variables(sm))

  # Agglomerated by phylum, it gives the coverage of the SingleM matrix
  phylum <- phyloseq::tax_glom(sm, "phylum", NArm = TRUE)
  otu <- as(phyloseq::otu_table(phylum), "matrix")
  names <- as.character(phyloseq::tax_table(phylum)[, "phylum"])
  wide <- read_out(dir, "matrices/singlem_phylum_coverage.tsv")
  sample <- "T1-BB-O-1"
  expect_equal(sort(unname(otu[, sample])),
               sort(wide[wide$taxon %in% names, sample]), tolerance = 1e-6)

  dir <- withr::local_tempdir()
  export_tables(extdata("marine_honour"), dir, excel = FALSE)
  mags <- readRDS(file.path(dir, "phyloseq", "mags.rds"))
  expect_equal(phyloseq::ntaxa(mags), 281)
  expect_equal(phyloseq::nsamples(mags), 48)
  expect_true(all(c("domain", "phylum", "species") %in% phyloseq::rank_names(mags)))
  expect_false(file.exists(file.path(dir, "phyloseq", "singlem.rds")))
})

test_that("export_tables uses the dereplication clusters and the CoverM abundance", {
  # Karen's case: --dereplicate without --mimag, each sample mapped to its own
  # assembly by the pipeline, and every sample mapped to the MAGs afterwards
  dir <- withr::local_tempdir()
  file.copy(extdata("example_run", "mag_quality.tsv"), dir)
  file.copy(extdata("example_run", "mapping_summary.tsv"), dir)
  quality <- read_mag_quality(file.path(dir, "mag_quality.tsv"))
  reps <- quality$bin[quality$quality == "medium"]
  others <- setdiff(quality$bin, reps)
  expect_length(reps, 3)
  utils::write.table(
    data.frame(representative = c(reps, reps[1], reps[1]),
               member = c(reps, others[1:2])),
    file.path(dir, "dereplication_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  samples <- c("minigut1", "minigut2")

  # The pipeline: each MAG only in one sample
  dir.create(file.path(dir, "04_binning"))
  own <- expand.grid(bin = quality$bin, sample = samples, stringsAsFactors = FALSE)
  own$mean_depth <- ifelse(own$sample == "minigut1", 5, 0)
  utils::write.table(own, file.path(dir, "04_binning", "mag_abundance.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  # CoverM: every MAG in both samples
  all <- expand.grid(bin = reps, sample = samples, stringsAsFactors = FALSE)
  all$mean_depth <- c(10, 4, 1, 2, 8, 3)
  all$covered_fraction <- 0.9
  all$reads_pct <- c(20, 8, 2, 4, 16, 6)
  dir.create(file.path(dir, "coverm"))
  utils::write.table(all, file.path(dir, "coverm", "mag_abundance_all_samples.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  run <- read_pipeline(dir)
  expect_identical(basename(attr(run, "files")$abundance), "mag_abundance_all_samples.tsv")
  expect_setequal(unique(run$abundance$bin), reps)
  expect_named(run$clusters, c("representative", "member"))

  out <- withr::local_tempdir()
  export_tables(run, out, excel = FALSE, phyloseq = FALSE)

  mags <- read_out(out, "mags.tsv")
  expect_setequal(mags$bin, reps)
  expect_identical(mags$n_cluster_members[match(reps, mags$bin)], c(3L, 1L, 1L))
  expect_equal(mags$reads_max_pct[mags$bin == reps[1]], 20)
  bins <- read_out(out, "bins.tsv")
  expect_identical(sum(bins$is_representative), 3L)

  reads <- read_out(out, "matrices/mag_reads_pct.tsv")
  expect_equal(dim(reads), c(3, 3))
  expect_true(file.exists(file.path(out, "matrices", "mag_covered_fraction.tsv")))
  sample_table <- read_out(out, "samples.tsv")
  expect_equal(sample_table$reads_in_mags_pct, c(30, 26))
  expect_true(file.exists(file.path(out, "long", "clusters.tsv")))

  readme <- readLines(file.path(out, "README.txt"))
  expect_false(any(grepl("one sample only", readme)))
  expect_true(any(grepl("every sample mapped to the dereplicated MAGs", readme)))
})
