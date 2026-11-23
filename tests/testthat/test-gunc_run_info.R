gunc_paths <- function() {
  list.files(system.file("extdata", "marine_honour", "gunc", package = "magreport"),
             full.names = TRUE)
}

test_that("read_gunc reads every report with its assembly and binner", {
  gunc <- read_gunc(gunc_paths())

  expect_equal(nrow(gunc), 1161)
  expect_setequal(unique(gunc$assembly), paste0("T", 0:7))
  expect_identical(unique(gunc$binner), "dastool")
  expect_type(gunc$pass, "logical")
  expect_false(anyNA(gunc$pass))
  expect_equal(sum(!gunc$pass), 157)
  expect_identical(gunc$bin, gunc$genome)
  expect_type(gunc$clade_separation_score, "double")
  expect_false("pass.GUNC" %in% names(gunc))
  # GUNC fails bins above its clade separation threshold (0.45); the score is
  # written with two decimals, so some failed bins show exactly 0.45
  expect_true(all(gunc$clade_separation_score[!gunc$pass] >= 0.45))
  expect_true(all(gunc$clade_separation_score[gunc$pass] <= 0.45))
})

test_that("read_gunc reads gzipped reports and names without a bin set", {
  src <- gunc_paths()[1]
  dir <- withr::local_tempdir()
  gz <- file.path(dir, "coassembly.gunc.tsv.gz")
  con <- gzfile(gz, "wb")
  writeLines(readLines(src), con)
  close(con)

  gunc <- read_gunc(gz)
  expect_equal(nrow(gunc), nrow(read_gunc(src)))
  expect_identical(unique(gunc$assembly), "coassembly")
  expect_true(all(is.na(gunc$binner)))
  expect_error(read_gunc(character()), "No GUNC report")
})

test_that("plot_gunc draws against CheckM2 or GUNC contamination", {
  skip_if_not_installed("ggplot2")
  gunc <- read_gunc(gunc_paths())
  mimag <- read_mimag_quality(system.file("extdata", "marine_honour", "mimag_quality.tsv",
                                          package = "magreport"))

  p <- plot_gunc(gunc)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "GUNC contamination portion")
  expect_match(p$labels$subtitle, "157 of 1161")

  # Only the 281 dereplicated MAGs have CheckM2 values in the MIMAG table
  p <- plot_gunc(gunc, mags = mimag)
  expect_identical(p$labels$x, "CheckM2 contamination (%)")
  expect_equal(nrow(p$data), 281)
  expect_match(p$labels$subtitle, "0 of 281")

  expect_error(plot_gunc(gunc, mags = data.frame(bin = "x", contamination = 1)), "No bin")
})

write_run_info <- function(env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".txt", .local_envir = env)
  writeLines(c(
    "mag-pipeline run information",
    "",
    "pipeline         : AntonioPuriel/mag-pipeline 2.0.0",
    "nextflow         : 26.04.6",
    "started          : 2026-11-11T09:17:54+01:00",
    "duration         : 1d 20h 40m 12s",
    "success          : true",
    "profile          : pyrene",
    "container engine : none (conda or local tools)",
    "",
    "command line:",
    "  nextflow run /path/mag-pipeline -profile pyrene --workflow both",
    "",
    "reference databases:",
    "  CheckM2 : /db/uniref100.KO.1.dmnd (last modified 2021-03-23)",
    "  GTDB-Tk : /db/gtdbtk (last modified 2026-05-31)",
    "  CAT     : not used",
    "  SingleM : /db/singlem (last modified unknown)",
    "",
    "parameters:",
    "  assembly_mode           = group",
    "  binners                 = concoct,metabat2",
    "  cat_db                  = null",
    "  catalog_chunk_size      = 1000000",
    "  catalog_min_id          = 0.95",
    "  refine                  = true",
    "  skip_fastp              = false"
  ), path)
  path
}

test_that("read_run_info reads the sections of run_info.txt", {
  info <- read_run_info(write_run_info())

  expect_identical(info$pipeline$name, "AntonioPuriel/mag-pipeline")
  expect_identical(info$pipeline$version, "2.0.0")
  expect_identical(info$pipeline$nextflow, "26.04.6")
  expect_true(info$run$success)
  expect_identical(info$run$container_engine, "none (conda or local tools)")
  expect_match(info$run$command_line, "--workflow both$")
  expect_identical(names(info$databases), c("checkm2", "gtdbtk", "singlem"))
  expect_identical(info$databases$gtdbtk$modified, "2026-05-31")
  expect_identical(info$databases$gtdbtk$path, "/db/gtdbtk")
  expect_null(info$databases$singlem$modified)
  expect_identical(info$params$catalog_chunk_size, 1e6)
  expect_identical(info$params$catalog_min_id, 0.95)
  expect_true(info$params$refine)
  expect_false(info$params$skip_fastp)
  expect_null(info$params$cat_db)
  expect_identical(info$params$binners, "concoct,metabat2")
})

test_that("describe_run gives the tables of the methods", {
  info <- read_run_info(write_run_info())
  methods <- describe_run(info)

  expect_identical(names(methods), c("run", "databases", "params"))
  expect_identical(methods$run$value[methods$run$item == "Version"], "2.0.0")
  expect_equal(nrow(methods$databases), 3)
  expect_true(is.na(methods$databases$modified[methods$databases$database == "singlem"]))
  # null parameters are left out; numbers are not written in scientific notation
  expect_false("cat_db" %in% methods$params$parameter)
  expect_identical(methods$params$value[methods$params$parameter == "catalog_chunk_size"], "1000000")

  expect_error(describe_run(list(a = 1)), "No run information")
})

test_that("read_pipeline keeps the run information and the GUNC reports", {
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "05_quality", "gunc"), recursive = TRUE)
  file.copy(gunc_paths(), file.path(dir, "05_quality", "gunc"))
  file.copy(system.file("extdata", "marine_honour", "mimag_quality.tsv", package = "magreport"),
            file.path(dir, "05_quality"))
  file.copy(write_run_info(), file.path(dir, "run_info.txt"))

  run <- read_pipeline(dir)
  expect_equal(nrow(run$gunc), 1161)
  expect_null(attr(run, "manifest"))
  expect_identical(attr(run, "run_info")$pipeline$version, "2.0.0")
  expect_identical(describe_run(run)$run$value[1], "AntonioPuriel/mag-pipeline")
})
