#' Plot how many MAGs were recovered per taxon
#'
#' Counts MAGs at the chosen taxonomic rank. MAGs that GTDB-Tk could not
#' classify at that rank are kept together in an "Unclassified" group, which is
#' usually the most informative part of the figure.
#'
#' @param mags A data frame returned by [read_mag_quality()].
#' @param rank Taxonomic rank: one of `"domain"`, `"phylum"`, `"class"`,
#'   `"order"`, `"family"`, `"genus"` or `"species"`.
#' @param top Maximum number of taxa to show. The rest are grouped into
#'   "Other". Use `NULL` to show them all.
#' @param fill_by Optional column used to colour the bars, for example
#'   `"quality"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "mag_quality.tsv",
#'                     package = "magreport")
#' mags <- read_mag_quality(path)
#' plot_taxonomy(mags, rank = "phylum")
#' plot_taxonomy(mags, rank = "genus", fill_by = "quality")
plot_taxonomy <- function(mags, rank = "phylum", top = 15, fill_by = NULL) {
  check_ggplot2()
  stopifnot(is.data.frame(mags))

  rank <- match.arg(
    rank,
    c("domain", "phylum", "class", "order", "family", "genus", "species")
  )

  if (!rank %in% names(mags)) {
    stop(
      "Column '", rank, "' not found: was the pipeline run with --gtdbtk_db?",
      call. = FALSE
    )
  }

  taxon <- as.character(mags[[rank]])
  taxon[is.na(taxon) | taxon == ""] <- "Unclassified"

  if (all(taxon == "Unclassified")) {
    stop("No MAG has a taxonomic assignment at rank '", rank, "'", call. = FALSE)
  }

  mags$taxon <- taxon

  if (!is.null(top)) {
    ranking <- sort(table(mags$taxon[mags$taxon != "Unclassified"]), decreasing = TRUE)
    keep <- names(ranking)[seq_len(min(top, length(ranking)))]
    mags$taxon[!mags$taxon %in% c(keep, "Unclassified")] <- "Other"
  }

  counts <- sort(table(mags$taxon[!mags$taxon %in% c("Other", "Unclassified")]))
  levels_ordered <- c("Unclassified", "Other", names(counts))
  levels_ordered <- levels_ordered[levels_ordered %in% unique(mags$taxon)]
  mags$taxon <- factor(mags$taxon, levels = levels_ordered)

  p <- if (is.null(fill_by)) {
    ggplot2::ggplot(mags, ggplot2::aes(x = .data$taxon)) +
      ggplot2::geom_bar(fill = "#2C7FB8", width = 0.7)
  } else {
    if (!fill_by %in% names(mags)) {
      stop("Column not found in mags: ", fill_by, call. = FALSE)
    }
    p <- ggplot2::ggplot(mags, ggplot2::aes(x = .data$taxon, fill = .data[[fill_by]])) +
      ggplot2::geom_bar(width = 0.7)
    if (fill_by == "quality") {
      p <- p + ggplot2::scale_fill_manual(
        name = "Quality", values = quality_colours(), drop = FALSE
      )
    }
    p
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "MAGs",
      title = paste0("MAGs per ", rank)
    ) +
    ggplot2::theme_bw()
}
