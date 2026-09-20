#' Plot which MAGs carry the genes of each process
#'
#' Crosses the KEGG orthologues found in every MAG with a set of marker genes
#' and draws a heatmap of MAGs against processes. It answers the question that
#' usually follows MAG recovery in a polluted site: which genomes can do what.
#'
#' @param kos A data frame returned by [read_ko_table()].
#' @param markers A data frame with the columns `pathway` and `ko`. Defaults to
#'   [hydrocarbon_kos()].
#' @param mags Optional data frame returned by [read_mag_quality()], used to
#'   label each MAG with its taxonomy instead of its bin name.
#' @param rank Taxonomic rank used for the labels when `mags` is given.
#' @param min_kos Minimum number of marker genes a MAG must carry to be shown.
#'   Use `0` to keep every MAG.
#'
#' @return A `ggplot` object. The fill is the number of marker genes of that
#'   process found in the MAG.
#' @export
#'
#' @examples
#' kos <- read_ko_table(system.file("extdata", "example_run", "ko_per_mag.tsv",
#'                                  package = "magreport"))
#' plot_pathway_presence(kos)
plot_pathway_presence <- function(kos,
                                  markers = hydrocarbon_kos(),
                                  mags = NULL,
                                  rank = "genus",
                                  min_kos = 1) {
  check_ggplot2()
  stopifnot(is.data.frame(kos), is.data.frame(markers))

  missing <- setdiff(c("pathway", "ko"), names(markers))
  if (length(missing) > 0) {
    stop("markers needs the columns pathway and ko", call. = FALSE)
  }

  hits <- merge(kos, markers, by = "ko")
  if (nrow(hits) == 0) {
    stop("None of the marker genes was found in these MAGs", call. = FALSE)
  }

  hits$n_genes <- if ("n_genes" %in% names(hits)) hits$n_genes else 1

  counts <- stats::aggregate(
    hits$n_genes,
    by = list(bin = hits$bin, binner = hits$binner, pathway = hits$pathway),
    FUN = sum
  )
  names(counts)[names(counts) == "x"] <- "n_genes"

  # Drop MAGs with too few marker genes to be informative
  per_bin <- stats::aggregate(counts$n_genes, by = list(bin = counts$bin), FUN = sum)
  keep <- per_bin$bin[per_bin$x >= min_kos]
  counts <- counts[counts$bin %in% keep, , drop = FALSE]

  if (nrow(counts) == 0) {
    stop("No MAG reaches min_kos = ", min_kos, call. = FALSE)
  }

  counts$label <- counts$bin
  if (!is.null(mags)) {
    if (!rank %in% names(mags)) {
      stop("Column '", rank, "' not found in mags", call. = FALSE)
    }
    taxon <- stats::setNames(as.character(mags[[rank]]), mags$bin)
    named <- taxon[counts$bin]
    named[is.na(named) | named == ""] <- NA
    counts$label <- ifelse(is.na(named), counts$bin, paste0(named, " (", counts$bin, ")"))
  }

  counts$pathway <- factor(counts$pathway, levels = unique(markers$pathway))

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = .data$pathway, y = .data$label, fill = .data$n_genes)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_genes), size = 3, colour = "grey15") +
    ggplot2::scale_fill_gradient(name = "Genes", low = "#DEEBF7", high = "#08519C") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Marker genes per MAG",
      caption = "Empty cells: no gene of that process was annotated"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}
