#' Read the GUNC chimera reports
#'
#' Reads the GUNC reports of mag-pipeline (`05_quality/gunc/<id>_<set>.gunc.tsv`,
#' written with `--gunc_db`), one per assembly and bin set. GUNC tells whether
#' a bin mixes genes of different lineages, which CheckM2 can miss when the
#' mixed genomes are incomplete. The pipeline runs it on every final bin,
#' before dereplication, so the reports also hold the bins that were left out
#' for being chimeric (`--gunc_exclude`).
#'
#' @param path One or more GUNC reports.
#'
#' @return A data frame with one row per bin: the GUNC columns, `bin` (a copy
#'   of `genome`, as in [read_mag_quality()]), `assembly` and `binner` (from the
#'   file name) and `pass` (logical, from `pass.GUNC`). A bin fails when its
#'   clade separation score (`clade_separation_score`) is above 0.45 at the
#'   taxonomic level where it is highest (`taxonomic_level`). GUNC writes the
#'   score with two decimals, so a failed bin can show exactly 0.45.
#' @export
#'
#' @examples
#' paths <- list.files(
#'   system.file("extdata", "marine_honour", "gunc", package = "magreport"),
#'   full.names = TRUE
#' )
#' gunc <- read_gunc(paths)
#' table(gunc$assembly, gunc$pass)
read_gunc <- function(path) {
  if (length(path) == 0) {
    stop("No GUNC report given", call. = FALSE)
  }

  one <- function(p) {
    x <- read_pipeline_table(
      p,
      required = c("genome", "taxonomic_level", "clade_separation_score", "pass.GUNC")
    )
    set <- sub("\\.gunc\\.tsv(\\.gz)?$", "", basename(p))
    x$assembly <- if (grepl("_", set)) sub("_[^_]+$", "", set) else set
    x$binner <- if (grepl("_", set)) sub(".*_", "", set) else NA_character_
    x
  }
  x <- do.call(rbind, lapply(path, one))

  numeric_cols <- c(
    "n_genes_called", "n_genes_mapped", "n_contigs",
    "proportion_genes_retained_in_major_clades", "genes_retained_index",
    "clade_separation_score", "contamination_portion",
    "n_effective_surplus_clades", "mean_hit_identity",
    "reference_representation_score"
  )
  for (column in intersect(numeric_cols, names(x))) {
    x[[column]] <- suppressWarnings(as.numeric(x[[column]]))
  }

  x$bin <- x$genome
  x$pass <- as.logical(x$pass.GUNC)
  x$pass.GUNC <- NULL
  rownames(x) <- NULL
  x
}

#' Plot GUNC chimerism against CheckM2 contamination
#'
#' Each point is a bin: the CheckM2 contamination on the x axis and the GUNC
#' clade separation score on the y axis, with the GUNC threshold (0.45). Bins
#' above the line are chimeric even when CheckM2 sees little contamination.
#' Without CheckM2 results, the x axis is the GUNC contamination portion.
#'
#' @param gunc A data frame returned by [read_gunc()].
#' @param mags Optional data frame returned by [read_mag_quality()] or
#'   [read_mimag_quality()], joined to `gunc` by `bin`. Its `quality` colours
#'   nothing here; it only provides the CheckM2 contamination.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' dir <- system.file("extdata", "marine_honour", package = "magreport")
#' gunc <- read_gunc(list.files(file.path(dir, "gunc"), full.names = TRUE))
#' plot_gunc(gunc)
#' plot_gunc(gunc, mags = read_mimag_quality(file.path(dir, "mimag_quality.tsv")))
plot_gunc <- function(gunc, mags = NULL) {
  check_ggplot2()
  stopifnot(is.data.frame(gunc))
  missing <- setdiff(c("bin", "clade_separation_score", "pass"), names(gunc))
  if (length(missing) > 0) {
    stop("Column(s) not found in gunc: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  x_label <- "GUNC contamination portion"
  gunc$x <- gunc$contamination_portion
  if (!is.null(mags) && all(c("bin", "contamination") %in% names(mags))) {
    gunc$x <- mags$contamination[match(gunc$bin, mags$bin)]
    gunc <- gunc[!is.na(gunc$x), , drop = FALSE]
    x_label <- "CheckM2 contamination (%)"
  }
  if (nrow(gunc) == 0) {
    stop("No bin of gunc is in mags", call. = FALSE)
  }

  gunc$status <- factor(ifelse(gunc$pass, "pass", "fail"), levels = c("pass", "fail"))
  col <- magreport_colours()
  n_fail <- sum(gunc$status == "fail", na.rm = TRUE)

  ggplot2::ggplot(gunc, ggplot2::aes(x = .data$x, y = .data$clade_separation_score)) +
    ggplot2::geom_hline(yintercept = 0.45, linetype = "dashed", colour = col[["sand_mid"]]) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status), size = 2.2, alpha = 0.75) +
    ggplot2::scale_colour_manual(
      name = "GUNC", values = c(pass = col[["olive"]], fail = col[["terracotta"]]), drop = FALSE
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Clade separation score",
      title = "Chimerism (GUNC)",
      subtitle = sprintf("%d of %d bins fail", n_fail, nrow(gunc)),
      caption = "Dashed line: GUNC threshold (0.45)"
    ) +
    theme_magreport()
}
