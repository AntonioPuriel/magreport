#' Read the KEGG orthologues found in each MAG
#'
#' Reads `07_function/ko_per_mag.tsv`, written by mag-pipeline when it is run
#' with `--eggnog_db`.
#'
#' @param path Path to `ko_per_mag.tsv`.
#'
#' @return A data frame with one row per MAG and KEGG orthologue.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "ko_per_mag.tsv",
#'                     package = "magreport")
#' kos <- read_ko_table(path)
#' head(kos)
read_ko_table <- function(path) {
  x <- read_pipeline_table(path, required = c("assembly", "binner", "bin", "ko"))

  if ("n_genes" %in% names(x)) {
    x$n_genes <- as.numeric(x$n_genes)
  }

  x
}

#' Read the annotation statistics of each MAG
#'
#' Reads `07_function/annotation_stats.tsv`: genes called by Prodigal, genes
#' annotated by eggNOG-mapper and distinct KEGG orthologues per MAG.
#'
#' @param path Path to `annotation_stats.tsv`.
#'
#' @return A data frame with one row per MAG and an extra column
#'   `annotated_pct`, the percentage of called genes that were annotated.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "annotation_stats.tsv",
#'                     package = "magreport")
#' read_annotation_stats(path)
read_annotation_stats <- function(path) {
  x <- read_pipeline_table(
    path,
    required = c("assembly", "binner", "bin", "n_genes", "n_annotated")
  )

  for (column in intersect(c("n_genes", "n_annotated", "n_kos"), names(x))) {
    x[[column]] <- as.numeric(x[[column]])
  }

  x$annotated_pct <- ifelse(x$n_genes > 0, 100 * x$n_annotated / x$n_genes, NA_real_)

  x
}
