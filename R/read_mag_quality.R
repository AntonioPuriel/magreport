#' Read the final MAG table
#'
#' Reads `05_quality/mag_quality.tsv`, the table produced by mag-pipeline with
#' one row per bin assessed by CheckM2, including the taxonomy assigned by
#' GTDB-Tk when it was run.
#'
#' @param path Path to `mag_quality.tsv`.
#'
#' @return A data frame with one row per MAG. `quality` is a factor with the
#'   levels returned by [quality_levels()]. Taxonomic ranks are `NA` when
#'   taxonomy was not assigned.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "mag_quality.tsv",
#'                     package = "magreport")
#' mags <- read_mag_quality(path)
#' head(mags)
read_mag_quality <- function(path) {
  x <- read_pipeline_table(
    path,
    required = c("assembly", "binner", "bin", "completeness", "contamination", "quality")
  )

  numeric_cols <- c(
    "n_contigs", "total_length", "n50", "gc_pct",
    "completeness", "contamination", "closest_ani"
  )
  for (column in intersect(numeric_cols, names(x))) {
    x[[column]] <- as.numeric(x[[column]])
  }

  x$quality <- factor(x$quality, levels = quality_levels())

  x
}

#' Count MAGs per quality category
#'
#' @param mags A data frame returned by [read_mag_quality()].
#' @param by Optional column names used to group the counts, for example
#'   `"binner"` or `c("assembly", "binner")`.
#'
#' @return A data frame with the number of MAGs in each quality category.
#'   Categories with no MAGs are kept with a count of zero.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "mag_quality.tsv",
#'                     package = "magreport")
#' mags <- read_mag_quality(path)
#' count_quality(mags)
#' count_quality(mags, by = "binner")
count_quality <- function(mags, by = NULL) {
  stopifnot(is.data.frame(mags), "quality" %in% names(mags))

  missing <- setdiff(by, names(mags))
  if (length(missing) > 0) {
    stop("Column(s) not found in mags: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  mags$quality <- factor(mags$quality, levels = quality_levels())

  if (is.null(by)) {
    counts <- as.data.frame(table(quality = mags$quality), stringsAsFactors = FALSE)
  } else {
    groups <- c(lapply(by, function(column) mags[[column]]), list(quality = mags$quality))
    names(groups) <- c(by, "quality")
    counts <- as.data.frame(table(groups), stringsAsFactors = FALSE)
  }

  names(counts)[names(counts) == "Freq"] <- "n"
  counts$quality <- factor(counts$quality, levels = quality_levels())
  counts <- counts[order(counts$quality), ]
  rownames(counts) <- NULL
  counts
}
