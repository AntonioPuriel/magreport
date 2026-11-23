# A run laid out as mag-pipeline publishes it, built from the example data,
# with a manifest that lists its tables
local_manifest_run <- function(schema_version = "1.0", extra_tables = list(),
                               env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  example <- function(file) system.file("extdata", "example_run", file, package = "magreport")

  layout <- list(
    "02_assembly" = "coassembly.contigs.stats.tsv",
    "03_mapping"  = c("mapping_summary.tsv", "coassembly.depth.txt"),
    "04_binning"  = "bin_summary.tsv",
    "05_quality"  = "mag_quality.tsv",
    "07_function" = c("ko_per_mag.tsv", "annotation_stats.tsv")
  )
  for (folder in names(layout)) {
    dir.create(file.path(dir, folder))
    file.copy(example(layout[[folder]]), file.path(dir, folder))
  }

  tables <- c(list(
    assembly_stats   = "02_assembly/coassembly.contigs.stats.tsv",
    mapping_summary  = "03_mapping/mapping_summary.tsv",
    depth            = "03_mapping/coassembly.depth.txt",
    bin_summary      = "04_binning/bin_summary.tsv",
    mag_quality      = "05_quality/mag_quality.tsv",
    ko_per_mag       = "07_function/ko_per_mag.tsv",
    annotation_stats = "07_function/annotation_stats.tsv"
  ), extra_tables)

  manifest <- list(
    schema_version = schema_version,
    pipeline = list(name = "AntonioPuriel/mag-pipeline", version = "2.1.0", commit = NULL),
    run = list(name = "test_run", success = TRUE, profile = "test"),
    databases = list(gtdbtk = list(path = "/db/gtdbtk", modified = "2026-05-02")),
    params = list(workflow = "mags", min_contig_length = 500, refine = FALSE),
    tables = lapply(tables, I)
  )
  jsonlite::write_json(manifest, file.path(dir, "manifest.json"),
                       auto_unbox = TRUE, null = "null", pretty = TRUE)
  dir
}

test_that("read_manifest reads every section and makes paths absolute", {
  dir <- local_manifest_run()
  manifest <- read_manifest(dir)

  expect_identical(manifest$schema_version, "1.0")
  expect_identical(manifest$pipeline$version, "2.1.0")
  expect_null(manifest$pipeline$commit)
  expect_true(manifest$run$success)
  expect_identical(manifest$params$min_contig_length, 500L)
  expect_identical(manifest$databases$gtdbtk$modified, "2026-05-02")
  expect_true(all(file.exists(unlist(manifest$tables))))
  expect_identical(basename(manifest$tables$mag_quality), "mag_quality.tsv")
  expect_identical(attr(manifest, "dir"), normalizePath(dir, winslash = "/"))

  # The file itself works as well as its directory
  expect_identical(read_manifest(file.path(dir, "manifest.json")), manifest)
})

test_that("read_manifest checks schema_version", {
  expect_error(read_manifest(local_manifest_run("2.0")), "schema_version 2.0")
  expect_warning(read_manifest(local_manifest_run("1.7")), "newer")
  expect_error(read_manifest(local_manifest_run("uno")), "Invalid schema_version")
  expect_error(read_manifest(withr::local_tempdir()), "Manifest not found")
})

test_that("read_pipeline reads the tables the manifest lists", {
  dir <- local_manifest_run()
  run <- read_pipeline(dir)
  flat <- read_pipeline(system.file("extdata", "example_run", package = "magreport"))

  expect_false(is.null(attr(run, "manifest")))
  expect_null(attr(flat, "manifest"))
  for (name in c("mags", "bins", "mapping", "depth", "assembly", "kos", "annotation")) {
    expect_equal(run[[name]], flat[[name]], info = name)
  }
  expect_null(run$mimag)
  expect_null(run$singlem)
  expect_match(attr(run, "files")$mags, "05_quality/mag_quality.tsv$")
})

test_that("read_pipeline only reads what the manifest lists", {
  dir <- local_manifest_run()
  # A table on disk that the manifest does not list is not read
  file.copy(
    system.file("extdata", "marine_honour", "mimag_quality.tsv", package = "magreport"),
    file.path(dir, "05_quality")
  )
  expect_null(read_pipeline(dir)$mimag)
})

test_that("read_pipeline warns about listed files that are missing", {
  dir <- local_manifest_run(extra_tables = list(gunc = "05_quality/gunc/x.gunc.tsv",
                                                mimag_quality = "05_quality/mimag_quality.tsv"))
  expect_warning(run <- read_pipeline(dir), "mimag_quality.tsv")
  expect_null(run$mimag)
  expect_false(is.null(run$mags))
})

test_that("read_pipeline reads the compressed tables of report/", {
  dir <- local_manifest_run()
  plain <- read_pipeline(dir)

  # report/ as mag-pipeline writes it: every .tsv and .txt gzipped, without
  # the per-contig tables, and its own manifest
  report <- file.path(dir, "report")
  manifest <- jsonlite::fromJSON(file.path(dir, "manifest.json"), simplifyDataFrame = FALSE)
  manifest$tables$depth <- NULL
  manifest$tables <- lapply(manifest$tables, function(files) {
    I(vapply(files, function(rel) {
      dest <- file.path(report, paste0(rel, ".gz"))
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      con <- gzfile(dest, "wb")
      writeLines(readLines(file.path(dir, rel)), con)
      close(con)
      paste0(rel, ".gz")
    }, character(1), USE.NAMES = FALSE))
  })
  manifest$bundle <- list(excluded = I("depth"))
  jsonlite::write_json(manifest, file.path(report, "manifest.json"),
                       auto_unbox = TRUE, null = "null", pretty = TRUE)

  run <- read_pipeline(report)
  expect_null(run$depth)
  expect_identical(attr(run, "manifest")$bundle$excluded, "depth")
  for (name in c("mags", "bins", "mapping", "assembly", "kos", "annotation")) {
    expect_equal(run[[name]], plain[[name]], info = name)
  }
})
