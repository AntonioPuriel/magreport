#' Read the MIMAG quality table
#'
#' Reads `05_quality/mimag_quality.tsv`, written by mag-pipeline (version 1.8.0
#' or later) when it is run with `--mimag`. It has one row per genome, with
#' CheckM2 completeness and contamination, the rRNA and tRNA genes found, the
#' MIMAG category and the GTDB-Tk classification.
#'
#' The columns `bin` (a copy of `genome`) and `quality` (the MIMAG category as
#' a factor) are added, and the GTDB-Tk classification is split into one
#' column per rank, so the result can be passed to the same functions as
#' [read_mag_quality()]: [count_quality()], [plot_mag_quality()] and
#' [plot_taxonomy()].
#'
#' @param path Path to `mimag_quality.tsv`.
#'
#' @return A data frame with one row per genome. `quality` is a factor with the
#'   levels returned by [quality_levels()]. Ranks without a name in GTDB are
#'   `NA`.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "marine_honour", "mimag_quality.tsv",
#'                     package = "magreport")
#' mags <- read_mimag_quality(path)
#' count_quality(mags)
#' table(mags$limiting)
read_mimag_quality <- function(path) {
  x <- read_pipeline_table(
    path,
    required = c("genome", "completeness", "contamination", "mimag_quality")
  )

  numeric_cols <- c(
    "n_contigs", "total_length", "n50", "gc_pct", "completeness", "contamination",
    "rrna_5S", "rrna_16S", "rrna_23S", "rrna_partial",
    "trna_types", "trna_total", "closest_ani"
  )
  for (column in intersect(numeric_cols, names(x))) {
    x[[column]] <- suppressWarnings(as.numeric(x[[column]]))
  }

  x$bin <- x$genome
  x$quality <- factor(x$mimag_quality, levels = quality_levels())

  if ("classification" %in% names(x)) {
    ranks <- split_gtdb_classification(x$classification)
    # The pipeline already writes the domain (from GTDB-Tk or the rRNA models)
    for (rank in setdiff(names(ranks), names(x))) {
      x[[rank]] <- ranks[[rank]]
    }
  }

  x
}

#' @keywords internal
split_gtdb_classification <- function(classification) {
  prefixes <- c(
    domain = "d", phylum = "p", class = "c", order = "o",
    family = "f", genus = "g", species = "s"
  )

  out <- lapply(prefixes, function(prefix) {
    value <- regmatches(
      classification,
      regexpr(paste0("(^|;)", prefix, "__[^;]*"), classification)
    )
    # regmatches drops non-matches: put them back as NA
    full <- rep(NA_character_, length(classification))
    hit <- grepl(paste0("(^|;)", prefix, "__"), classification)
    full[hit] <- sub(paste0("^;?", prefix, "__"), "", value)
    full[!is.na(full) & full == ""] <- NA_character_
    full
  })

  as.data.frame(out, stringsAsFactors = FALSE)
}
