#' Plot which hydrocarbon degradation markers each MAG carries
#'
#' Draws a heatmap of MAGs against the markers found by CANT-HYD or HADEG,
#' with the number of genes in each cell. MAGs are ordered by the number of
#' different markers they carry, so the most versatile degraders come first.
#'
#' @param genes A data frame returned by [read_hydrocarbon_genes()].
#' @param database `"CANT-HYD"` or `"HADEG"`.
#' @param mags Optional data frame returned by [read_mimag_quality()] or
#'   [read_mag_quality()], used to label each MAG with its taxonomy.
#' @param rank Taxonomic rank used for the labels when `mags` is given. MAGs
#'   without a name at that rank take the closest named rank above.
#' @param top Maximum number of MAGs to show, or `NULL` for all of them. HADEG
#'   finds markers in most MAGs, so a limit keeps the figure readable.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' dir <- system.file("extdata", "marine_honour", package = "magreport")
#' genes <- read_hydrocarbon_genes(file.path(dir, "hydrocarbon_per_mag.tsv"))
#' mags <- read_mimag_quality(file.path(dir, "mimag_quality.tsv"))
#' plot_hydrocarbon_genes(genes, mags = mags)
plot_hydrocarbon_genes <- function(genes,
                                   database = c("CANT-HYD", "HADEG"),
                                   mags = NULL,
                                   rank = "genus",
                                   top = 40) {
  check_ggplot2()
  stopifnot(is.data.frame(genes))
  database <- match.arg(database)

  hits <- genes[genes$database == database, , drop = FALSE]
  if (nrow(hits) == 0) {
    stop("No ", database, " markers in genes", call. = FALSE)
  }

  counts <- stats::aggregate(n_genes ~ bin + marker, data = hits, FUN = sum)

  richness <- sort(table(counts$bin), decreasing = TRUE)
  keep <- names(richness)
  if (!is.null(top)) {
    keep <- utils::head(keep, top)
  }
  counts <- counts[counts$bin %in% keep, , drop = FALSE]

  counts$label <- counts$bin
  if (!is.null(mags)) {
    ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
    if (!rank %in% ranks || !rank %in% names(mags)) {
      stop("Column '", rank, "' not found in mags", call. = FALSE)
    }
    upper <- rev(ranks[seq_len(match(rank, ranks))])
    upper <- upper[upper %in% names(mags)]
    name <- rep(NA_character_, nrow(mags))
    for (r in upper) {
      fill <- is.na(name) & !is.na(mags[[r]]) & mags[[r]] != ""
      name[fill] <- if (r == rank) mags[[r]][fill] else paste0(mags[[r]][fill], " (", substr(r, 1, 1), ")")
    }
    taxon <- stats::setNames(name, mags$bin)[counts$bin]
    counts$label <- ifelse(is.na(taxon), counts$bin, paste0(taxon, "  ", counts$bin))
  }

  label_order <- unique(counts$label[order(match(counts$bin, keep))])
  counts$label <- factor(counts$label, levels = rev(label_order))

  marker_order <- names(sort(table(counts$marker), decreasing = TRUE))
  counts$marker <- factor(counts$marker, levels = marker_order)

  # A single MAG with many copies of one marker would wash out the rest of the
  # scale: colours stop at the 95th percentile, the numbers show the real count
  cap <- max(3, ceiling(stats::quantile(counts$n_genes, 0.95, names = FALSE)))
  counts$fill <- pmin(counts$n_genes, cap)

  # Light text on the dark end of the scale, dark text elsewhere
  dark_tile <- counts$fill >= 0.55 * cap
  counts$text_colour <- ifelse(dark_tile, "white", magreport_colours()[["olive_dark"]])

  ggplot2::ggplot(counts, ggplot2::aes(x = .data$marker, y = .data$label, fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_genes, colour = .data$text_colour), size = 2.6) +
    ggplot2::scale_colour_identity() +
    scale_fill_magreport_c(
      name = "Genes", limits = c(0, cap),
      breaks = c(0, cap), labels = c("0", paste0(cap, "+"))
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = paste(database, "markers per MAG"),
      caption = if (!is.null(top) && length(richness) > top) {
        paste0("The ", top, " MAGs with most markers, out of ", length(richness))
      } else {
        NULL
      }
    ) +
    theme_magreport() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}
