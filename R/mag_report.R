#' Build an HTML report of a mag-pipeline run
#'
#' Reads the tables found in a results directory and renders a self-contained
#' HTML report with the mapping, assembly, binning and MAG quality results.
#' Sections without data, such as taxonomy when GTDB-Tk was not run, are left
#' out instead of failing.
#'
#' @param results_dir Directory holding the output of a mag-pipeline run.
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
  if (is.null(run$mags) && is.null(run$bins)) {
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
    params = list(run = run, title = title, results_dir = normalizePath(results_dir)),
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )

  invisible(output_file)
}
