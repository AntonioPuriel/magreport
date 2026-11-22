#' Plot the composition of the community in every sample
#'
#' Draws stacked bars with the relative abundance of the main taxa of a rank in
#' every sample. It works with the profile of the gene catalog (CAT, from
#' [read_taxonomy_profile()]) and with the profile of the reads (SingleM, from
#' [read_singlem()]), so the two can be drawn side by side.
#'
#' @param profile A data frame returned by [read_taxonomy_profile()] or
#'   [read_singlem()].
#' @param rank Taxonomic rank: one of `"domain"`, `"phylum"`, `"class"`,
#'   `"order"`, `"family"`, `"genus"` or `"species"`.
#' @param top Number of taxa shown, the most abundant on average across
#'   samples. The rest are grouped into "Other". `"Unclassified"` is always
#'   shown apart.
#' @param samples Optional data frame with one row per sample: the column
#'   `sample` and any design columns, for example from [parse_sample_names()].
#'   Samples are drawn in the order of its rows.
#' @param facet_by Optional column of `samples` used to split the samples into
#'   panels, for example a site or a condition.
#' @param strip_prefix Remove the GTDB rank prefix (`p__`) from the names.
#' @param title Title of the plot. By default it names the rank and the
#'   source of the profile.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_contigs", "profiles",
#'                     "taxonomy_per_sample.tsv.gz", package = "magreport")
#' cat_profile <- read_taxonomy_profile(path)
#' plot_community(cat_profile, rank = "phylum")
#'
#' samples <- parse_sample_names(
#'   cat_profile$sample,
#'   into = c("time", "site", "condition", "replicate")
#' )
#' plot_community(cat_profile, rank = "genus", samples = samples, facet_by = "condition")
plot_community <- function(profile,
                           rank = "phylum",
                           top = 10,
                           samples = NULL,
                           facet_by = NULL,
                           strip_prefix = TRUE,
                           title = NULL) {
  check_ggplot2()
  check_profile(profile)
  rank <- match.arg(rank, taxonomic_ranks())
  stopifnot(is.numeric(top), length(top) == 1, top >= 1)

  x <- profile[profile$rank == rank, c("taxon", "sample", "rel_abundance"), drop = FALSE]
  if (nrow(x) == 0) {
    stop("No taxa at rank '", rank, "' in profile", call. = FALSE)
  }

  n_samples <- length(unique(x$sample))
  classified <- x[x$taxon != "Unclassified", , drop = FALSE]
  means <- sort(tapply(classified$rel_abundance, classified$taxon, sum) / n_samples,
                decreasing = TRUE)
  keep <- names(means)[seq_len(min(top, length(means)))]

  x$taxon[!x$taxon %in% c(keep, "Unclassified")] <- "Other"
  x <- stats::aggregate(rel_abundance ~ taxon + sample, data = x, FUN = sum)

  levels_ordered <- c(keep, "Other", "Unclassified")
  levels_ordered <- levels_ordered[levels_ordered %in% x$taxon]
  colours <- c(
    stats::setNames(magreport_palette(length(keep)), keep),
    Other = "#B9AD95",
    Unclassified = "#E6DFD3"
  )[levels_ordered]
  shown <- if (strip_prefix) strip_rank_prefix(levels_ordered) else levels_ordered
  x$taxon <- factor(x$taxon, levels = levels_ordered)

  x <- attach_samples(x, samples, facet_by)

  source <- attr(profile, "source")
  if (is.null(title)) {
    title <- paste0("Community composition: ", rank,
                    if (!is.null(source)) paste0(" (", source, ")") else "")
  }

  p <- ggplot2::ggplot(
    x,
    ggplot2::aes(x = .data$sample, y = .data$rel_abundance, fill = .data$taxon)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_fill_manual(name = NULL, values = colours, labels = shown, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(
      x = NULL,
      y = "Relative abundance (%)",
      title = title,
      caption = profile_caption(source)
    ) +
    theme_magreport() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = ggplot2::element_blank()
    )

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_grid(
      cols = ggplot2::vars(.data[[facet_by]]),
      scales = "free_x",
      space = "free_x"
    )
  }
  p
}

#' Compare the community seen by CAT and by SingleM
#'
#' Plots the mean relative abundance of every taxon in the profile of the gene
#' catalog against that in the profile of the reads. Both are first rescaled to
#' their classified part, so that the large unclassified fraction of the
#' contigs does not shift every point. Taxa on the diagonal are seen alike by
#' both methods; those far below it are missed by the assemblies.
#'
#' @param x,y Data frames returned by [read_taxonomy_profile()] or
#'   [read_singlem()], usually one of each.
#' @param rank Taxonomic rank compared.
#' @param labels Names of the two profiles, for the axes. By default their
#'   `"source"` attribute.
#' @param label_top Number of taxa named on the plot, the most abundant ones.
#' @param strip_prefix Remove the GTDB rank prefix from the names.
#'
#' @return A `ggplot` object. Its `data` element holds the mean abundance of
#'   every taxon in both profiles, and the subtitle the Spearman correlation.
#' @export
#'
#' @examples
#' dir <- system.file("extdata", "example_contigs", package = "magreport")
#' cat_profile <- read_taxonomy_profile(
#'   file.path(dir, "profiles", "taxonomy_per_sample.tsv.gz")
#' )
#' singlem <- read_singlem(list.files(file.path(dir, "singlem"), full.names = TRUE))
#' plot_profile_comparison(cat_profile, singlem, rank = "phylum")
plot_profile_comparison <- function(x,
                                    y,
                                    rank = "phylum",
                                    labels = NULL,
                                    label_top = 8,
                                    strip_prefix = TRUE) {
  check_ggplot2()
  check_profile(x)
  check_profile(y)
  rank <- match.arg(rank, taxonomic_ranks())

  if (is.null(labels)) {
    labels <- c(attr(x, "source") %||% "Profile x", attr(y, "source") %||% "Profile y")
  }
  stopifnot(length(labels) == 2)

  mean_classified <- function(p) {
    p <- p[p$rank == rank & p$taxon != "Unclassified", c("taxon", "sample", "rel_abundance")]
    if (nrow(p) == 0) {
      stop("No classified taxa at rank '", rank, "'", call. = FALSE)
    }
    p$share <- share_within(p$rel_abundance, p$sample)
    n <- length(unique(p$sample))
    sums <- tapply(p$share, p$taxon, sum) / n
    data.frame(taxon = names(sums), mean = as.numeric(sums), stringsAsFactors = FALSE)
  }

  shared <- intersect(unique(x$sample), unique(y$sample))
  if (length(shared) == 0) {
    stop("The two profiles have no sample in common", call. = FALSE)
  }
  a <- mean_classified(x[x$sample %in% shared, ])
  b <- mean_classified(y[y$sample %in% shared, ])
  both <- merge(a, b, by = "taxon", all = TRUE, suffixes = c("_x", "_y"))
  both$mean_x[is.na(both$mean_x)] <- 0
  both$mean_y[is.na(both$mean_y)] <- 0

  rho <- suppressWarnings(stats::cor(both$mean_x, both$mean_y, method = "spearman"))
  ranking <- order(-(both$mean_x + both$mean_y))
  both$label <- ""
  top <- ranking[seq_len(min(label_top, nrow(both)))]
  both$label[top] <- if (strip_prefix) strip_rank_prefix(both$taxon[top]) else both$taxon[top]

  limit <- max(c(both$mean_x, both$mean_y))
  col <- magreport_colours()

  ggplot2::ggplot(both, ggplot2::aes(x = .data$mean_x, y = .data$mean_y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = col[["taupe"]]) +
    ggplot2::geom_point(colour = col[["olive"]], size = 2.2, alpha = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      size = 3, hjust = -0.12, vjust = 0.5, colour = col[["olive_dark"]], check_overlap = TRUE
    ) +
    # room on the right for the names of the most abundant taxa
    ggplot2::scale_x_sqrt(limits = c(0, limit * 1.35)) +
    ggplot2::scale_y_sqrt(limits = c(0, limit * 1.05)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = paste0(labels[1], ": mean abundance (%)"),
      y = paste0(labels[2], ": mean abundance (%)"),
      title = paste0("Two views of the community: ", rank),
      subtitle = sprintf("%d taxa, %d samples. Spearman rho = %.2f",
                         nrow(both), length(shared), rho),
      caption = paste(
        "Abundance among the classified part of each profile, averaged over samples.",
        "Square-root axes. Dashed line: same abundance in both.",
        sep = "\n"
      )
    ) +
    theme_magreport()
}

#' Plot how much of each sample the annotations explain
#'
#' Draws, for every annotation of the gene catalog, the share of the TPM of
#' each sample carried by annotated genes: genes with a KEGG orthologue, a
#' hydrocarbon marker, or on a contig that CAT classified.
#'
#' @param stats A data frame returned by [read_profile_stats()].
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' stats <- read_profile_stats(system.file("extdata", "example_contigs", "profiles",
#'                                         "profile_stats.tsv", package = "magreport"))
#' plot_annotation_fraction(stats)
plot_annotation_fraction <- function(stats) {
  check_ggplot2()
  stopifnot(is.data.frame(stats))
  columns <- grep("^fraction_", names(stats), value = TRUE)
  if (!"sample" %in% names(stats) || length(columns) == 0) {
    stop("stats must come from read_profile_stats()", call. = FALSE)
  }

  names_known <- c(
    fraction_eggnog_ko = "KEGG orthologue (eggNOG)",
    fraction_canthyd   = "CANT-HYD marker",
    fraction_hadeg     = "HADEG gene",
    fraction_cat_domain = "Taxonomy (CAT, domain)"
  )
  shown <- ifelse(columns %in% names(names_known), names_known[columns],
                  sub("^fraction_", "", columns))

  long <- do.call(rbind, lapply(seq_along(columns), function(i) {
    data.frame(sample = stats$sample, source = shown[i],
               percent = 100 * stats[[columns[i]]], stringsAsFactors = FALSE)
  }))
  long$source <- factor(long$source, levels = rev(shown))
  col <- magreport_colours()

  ggplot2::ggplot(long, ggplot2::aes(x = .data$source, y = .data$percent)) +
    ggplot2::geom_boxplot(fill = col[["sand"]], colour = col[["taupe"]],
                          outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(colour = col[["olive"]], width = 0.15, height = 0,
                         size = 1.6, alpha = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::expand_limits(y = c(0, 100)) +
    ggplot2::labs(
      x = NULL,
      y = "Share of the sample's TPM (%)",
      title = "Annotated part of the gene catalog",
      caption = "One point per sample"
    ) +
    theme_magreport()
}

#' Plot the abundance of genes or markers in every sample
#'
#' Draws a heatmap of the abundance (TPM) of KEGG orthologues or hydrocarbon
#' markers of the gene catalog across samples.
#'
#' @param profile A data frame returned by [read_feature_profile()].
#' @param features Features to show: a character vector, or a data frame
#'   with the columns `ko` (or `feature`) and `pathway`, such as
#'   [hydrocarbon_kos()], to group them in panels by pathway. By default the
#'   `top` most abundant features.
#' @param top Number of features shown when `features` is not given.
#' @param samples,facet_by Optional sample table and column to order and split
#'   the samples, as in [plot_community()].
#'
#' @return A `ggplot` object. The fill is log10(TPM + 1); empty cells were not
#'   detected in that sample.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_contigs", "profiles",
#'                     "canthyd_per_sample.tsv.gz", package = "magreport")
#' markers <- read_feature_profile(path)
#' plot_feature_heatmap(markers)
plot_feature_heatmap <- function(profile,
                                 features = NULL,
                                 top = 25,
                                 samples = NULL,
                                 facet_by = NULL) {
  check_ggplot2()
  stopifnot(is.data.frame(profile))
  if (!all(c("feature", "sample", "tpm") %in% names(profile))) {
    stop("profile must come from read_feature_profile()", call. = FALSE)
  }

  pathway <- NULL
  if (is.data.frame(features)) {
    id <- intersect(c("feature", "ko", "marker"), names(features))[1]
    if (is.na(id) || !"pathway" %in% names(features)) {
      stop("features needs the columns ko (or feature) and pathway", call. = FALSE)
    }
    pathway <- stats::setNames(as.character(features$pathway), features[[id]])
    pathway <- pathway[!duplicated(names(pathway))]
    features <- names(pathway)
  }

  x <- stats::aggregate(tpm ~ feature + sample, data = profile, FUN = sum)
  if (is.null(features)) {
    means <- sort(tapply(x$tpm, x$feature, sum), decreasing = TRUE)
    features <- names(means)[seq_len(min(top, length(means)))]
  }
  found <- intersect(features, x$feature)
  if (length(found) == 0) {
    stop("None of the features was found in profile", call. = FALSE)
  }

  grid <- expand.grid(feature = found, sample = unique(profile$sample),
                      stringsAsFactors = FALSE)
  x <- merge(grid, x, by = c("feature", "sample"), all.x = TRUE)
  x$fill <- ifelse(is.na(x$tpm) | x$tpm <= 0, NA_real_, log10(x$tpm + 1))
  x$feature <- factor(x$feature, levels = rev(found))
  if (!is.null(pathway)) {
    x$pathway <- factor(pathway[as.character(x$feature)], levels = unique(pathway[found]))
  }

  x <- attach_samples(x, samples, facet_by)

  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$sample, y = .data$feature, fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_magreport_c(name = "log10(TPM + 1)", na.value = "#FBF8F2") +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Abundance in the gene catalog",
      caption = "Empty cells: not detected in that sample"
    ) +
    theme_magreport() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )

  rows <- if (!is.null(pathway)) ggplot2::vars(.data$pathway) else NULL
  cols <- if (!is.null(facet_by)) ggplot2::vars(.data[[facet_by]]) else NULL
  if (!is.null(rows) || !is.null(cols)) {
    p <- p + ggplot2::facet_grid(rows = rows, cols = cols, scales = "free", space = "free") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
  }
  p
}

# Shared helpers -------------------------------------------------------------

check_profile <- function(profile) {
  stopifnot(is.data.frame(profile))
  if (!all(c("rank", "taxon", "sample", "rel_abundance") %in% names(profile))) {
    stop("profile must come from read_taxonomy_profile() or read_singlem()", call. = FALSE)
  }
  invisible(TRUE)
}

strip_rank_prefix <- function(taxon) {
  sub("^[dpcofgs]__", "", taxon)
}

profile_caption <- function(source) {
  if (identical(source, "CAT")) {
    "Share of the TPM of the gene catalog; genes take the CAT classification of their contig"
  } else if (identical(source, "SingleM")) {
    "Share of the coverage of single-copy marker genes in the reads (SingleM)"
  } else {
    NULL
  }
}

# Orders the samples as in `samples` and adds the column used for panels
attach_samples <- function(x, samples, facet_by) {
  if (is.null(samples)) {
    if (!is.null(facet_by)) {
      stop("facet_by needs the samples table", call. = FALSE)
    }
    x$sample <- factor(x$sample, levels = sort(unique(as.character(x$sample))))
    return(x)
  }
  stopifnot(is.data.frame(samples))
  if (!"sample" %in% names(samples)) {
    stop("samples needs the column sample", call. = FALSE)
  }
  if (!is.null(facet_by) && !facet_by %in% names(samples)) {
    stop("Column not found in samples: ", facet_by, call. = FALSE)
  }
  absent <- setdiff(unique(as.character(x$sample)), samples$sample)
  if (length(absent) > 0) {
    stop(length(absent), " sample(s) are missing from samples: ",
         paste(utils::head(absent, 3), collapse = ", "), call. = FALSE)
  }
  order <- unique(as.character(samples$sample))
  if (!is.null(facet_by)) {
    x[[facet_by]] <- samples[[facet_by]][match(as.character(x$sample), samples$sample)]
    x[[facet_by]] <- factor(x[[facet_by]], levels = unique(samples[[facet_by]]))
  }
  x$sample <- factor(x$sample, levels = order[order %in% x$sample])
  x
}

`%||%` <- function(a, b) if (is.null(a)) b else a
