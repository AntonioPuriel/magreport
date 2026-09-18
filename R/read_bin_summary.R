#' Read the per-bin statistics
#'
#' Reads `04_binning/bin_summary.tsv`, with one row per bin produced by each
#' binner, including bins that CheckM2 did not assess.
#'
#' @param path Path to `bin_summary.tsv`.
#'
#' @return A data frame with one row per bin.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "bin_summary.tsv",
#'                     package = "magreport")
#' bins <- read_bin_summary(path)
#' head(bins)
read_bin_summary <- function(path) {
  x <- read_pipeline_table(
    path,
    required = c("assembly", "binner", "bin", "n_contigs", "total_length")
  )

  for (column in intersect(c("n_contigs", "total_length", "n50", "gc_pct"), names(x))) {
    x[[column]] <- as.numeric(x[[column]])
  }

  x
}
