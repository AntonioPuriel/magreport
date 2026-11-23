#' Share of the community represented by the MAGs
#'
#' For each sample and rank, the share of a community profile (SingleM, from
#' the reads, or CAT, from the gene catalog) that belongs to taxa with at
#' least one MAG. A genus at 30 % in SingleM that has a MAG counts its 30 %;
#' one without MAG counts nothing. It measures how much of the community the
#' MAGs cover taxonomically, not how many reads map to them.
#'
#' Taxa are matched by name, so the profile and the MAGs should use the same
#' taxonomy: GTDB for both SingleM and GTDB-Tk. Names that changed between the
#' GTDB release of the SingleM metapackage and that of GTDB-Tk are not
#' matched, which underestimates the share, more at the genus than at the
#' phylum level.
#'
#' @param profile A data frame returned by [read_singlem()] or
#'   [read_taxonomy_profile()].
#' @param mags A data frame returned by [read_mimag_quality()] or
#'   [read_mag_quality()], with one column per rank.
#' @param ranks Ranks to compute.
#' @param min_quality Lowest quality of the MAGs that count: `"medium"` (the
#'   default) leaves low-quality bins out, `"low"` counts every MAG.
#'
#' @return A data frame with one row per sample and rank: `sample`, `rank`,
#'   `n_taxa_with_mag` (taxa of the profile with a MAG in that sample),
#'   `classified` (% of the profile classified at that rank) and `recovered`
#'   (% of the profile in taxa with a MAG). `recovered_of_classified` is
#'   `recovered` as a share of `classified`.
#' @export
#'
#' @examples
#' singlem <- read_singlem(list.files(
#'   system.file("extdata", "example_contigs", "singlem", package = "magreport"),
#'   full.names = TRUE
#' ))
#' mags <- data.frame(
#'   bin = c("a", "b"),
#'   quality = c("high", "medium"),
#'   phylum = c("Bacteroidota", "Pseudomonadota"),
#'   genus = c("Polaribacter", "Alcanivorax")
#' )
#' recovery <- mag_recovery(singlem, mags, ranks = c("phylum", "genus"))
#' head(recovery)
mag_recovery <- function(profile, mags,
                         ranks = c("phylum", "class", "order", "family", "genus"),
                         min_quality = c("medium", "high", "low")) {
  check_profile(profile)
  stopifnot(is.data.frame(mags))
  min_quality <- match.arg(min_quality)

  ranks <- intersect(ranks, names(mags))
  if (length(ranks) == 0) {
    stop("mags has no taxonomic rank columns (phylum, class, ...): run GTDB-Tk", call. = FALSE)
  }
  if ("quality" %in% names(mags)) {
    keep <- quality_levels()[seq_len(match(min_quality, quality_levels()))]
    mags <- mags[!is.na(mags$quality) & as.character(mags$quality) %in% keep, , drop = FALSE]
  }

  one_rank <- function(rank) {
    x <- profile[profile$rank == rank, c("taxon", "sample", "rel_abundance"), drop = FALSE]
    if (nrow(x) == 0) {
      return(NULL)
    }
    with_mag <- unique(stats::na.omit(strip_rank_prefix(mags[[rank]])))
    x$name <- strip_rank_prefix(x$taxon)
    x$classified <- x$taxon != "Unclassified" & x$name != ""
    x$has_mag <- x$classified & x$name %in% with_mag

    samples <- sort(unique(x$sample))
    out <- data.frame(
      sample = samples,
      rank = rank,
      n_taxa_with_mag = as.integer(tapply(x$has_mag, factor(x$sample, samples), sum)),
      classified = as.numeric(tapply(x$rel_abundance * x$classified, factor(x$sample, samples), sum)),
      recovered = as.numeric(tapply(x$rel_abundance * x$has_mag, factor(x$sample, samples), sum)),
      stringsAsFactors = FALSE
    )
    out$recovered_of_classified <- ifelse(out$classified > 0, 100 * out$recovered / out$classified, NA)
    out
  }

  out <- do.call(rbind, lapply(ranks, one_rank))
  if (is.null(out)) {
    stop("The profile has none of the ranks ", paste(ranks, collapse = ", "), call. = FALSE)
  }
  out$rank <- factor(out$rank, levels = ranks)
  rownames(out) <- NULL
  attr(out, "source") <- attr(profile, "source")
  out
}

#' Plot the share of the community represented by the MAGs
#'
#' One point per sample at each rank: the share of the community profile in
#' taxa with a MAG ([mag_recovery()]), with a box over the samples. The grey
#' marks show the share classified at that rank, the most the MAGs could
#' reach.
#'
#' @param recovery A data frame returned by [mag_recovery()].
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' singlem <- read_singlem(list.files(
#'   system.file("extdata", "example_contigs", "singlem", package = "magreport"),
#'   full.names = TRUE
#' ))
#' mags <- data.frame(
#'   bin = c("a", "b"),
#'   phylum = c("Bacteroidota", "Pseudomonadota"),
#'   genus = c("Polaribacter", "Alcanivorax")
#' )
#' plot_mag_recovery(mag_recovery(singlem, mags, ranks = c("phylum", "genus")))
plot_mag_recovery <- function(recovery) {
  check_ggplot2()
  stopifnot(is.data.frame(recovery))
  missing <- setdiff(c("sample", "rank", "classified", "recovered"), names(recovery))
  if (length(missing) > 0) {
    stop("recovery must come from mag_recovery()", call. = FALSE)
  }
  col <- magreport_colours()
  source <- attr(recovery, "source")

  ggplot2::ggplot(recovery, ggplot2::aes(x = .data$rank)) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$classified),
      shape = 95, size = 8, colour = col[["sand_mid"]]
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(y = .data$recovered),
      width = 0.5, outlier.shape = NA, fill = col[["olive"]], alpha = 0.25,
      colour = col[["olive_dark"]]
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(y = .data$recovered),
      width = 0.12, height = 0, size = 1.8, alpha = 0.7, colour = col[["olive_dark"]]
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(
      x = NULL,
      y = "Share of the community (%)",
      title = "Community represented by the MAGs",
      subtitle = paste0(
        "Share of the ", if (is.null(source)) "profile" else source,
        " profile in taxa with at least one MAG; one point per sample"
      ),
      caption = "Grey marks: share classified at that rank"
    ) +
    theme_magreport()
}
