#' Build an HTML report of a mag-pipeline run
#'
#' Reads the tables found in a results directory and renders a self-contained
#' HTML report with the mapping, assembly, binning and MAG quality results,
#' and the profile of the community from the gene catalog and the reads (CAT,
#' eggNOG, CANT-HYD and SingleM; mag-pipeline 2.0.0 or later), the GUNC chimera
#' check and the methods of the run (versions, databases and parameters, from
#' `manifest.json` or `run_info.txt`). Sections without
#' data, such as taxonomy when GTDB-Tk was not run, are left out instead of
#' failing, so a run with `--workflow contigs` and no MAGs works too.
#'
#' @param results_dir Directory holding the output of a mag-pipeline run, or
#'   several, named by study, to compare them (see [read_pipeline()]). With
#'   several studies the report adds a table of the studies, splits the
#'   quality plots by study and leaves out the per-sample composition plots.
#' @param output_file Path of the HTML file to write.
#' @param title Title shown at the top of the report.
#' @param quiet Passed to [rmarkdown::render()]; `TRUE` hides the progress output.
#'
#' @return The path of the report, invisibly.
#' @export
#'
#' @examples
#' \donttest{
#' results_dir <- system.file("extdata", "example_run", package = "magreport")
#' if (rmarkdown::pandoc_available()) {
#'   mag_report(results_dir, file.path(tempdir(), "mag_report.html"))
#' }
#' }
mag_report <- function(results_dir,
                       output_file = "mag_report.html",
                       title = "MAG recovery report",
                       quiet = TRUE) {
  for (pkg in c("rmarkdown", "knitr", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is needed to build the report.", call. = FALSE)
    }
  }
  if (!rmarkdown::pandoc_available()) {
    stop("Pandoc is needed to build the report. It ships with RStudio.", call. = FALSE)
  }

  run <- read_pipeline(results_dir)
  found <- c("mags", "bins", "mimag", "catalog", "profile_stats", "contig_taxonomy", "singlem")
  if (all(vapply(run[found], is.null, logical(1)))) {
    stop("No mag-pipeline tables found in ", results_dir, call. = FALSE)
  }

  template <- system.file("rmd", "mag_report.Rmd", package = "magreport")
  if (template == "") {
    stop("Report template not found in the installed package", call. = FALSE)
  }

  output_file <- normalizePath(output_file, mustWork = FALSE)

  rmarkdown::render(
    input = template,
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    intermediates_dir = tempdir(),
    knit_root_dir = tempdir(),
    params = list(run = run, title = title,
                  results_dir = paste(normalizePath(results_dir), collapse = ", ")),
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )

  invisible(output_file)
}
