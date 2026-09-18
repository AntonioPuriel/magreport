#' Read the read mapping summary
#'
#' Reads `03_mapping/mapping_summary.tsv`, with the number of read pairs and
#' the overall Bowtie2 alignment rate of each sample.
#'
#' @param path Path to `mapping_summary.tsv`.
#'
#' @return A data frame with one row per sample.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "mapping_summary.tsv",
#'                     package = "magreport")
#' read_mapping_summary(path)
read_mapping_summary <- function(path) {
  x <- read_pipeline_table(path, required = c("sample", "read_pairs"))

  for (column in intersect(c("read_pairs", "overall_alignment_rate_pct"), names(x))) {
    x[[column]] <- as.numeric(x[[column]])
  }

  x
}
