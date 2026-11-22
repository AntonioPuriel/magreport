#' Plot the abundance of groups of MAGs, genes or taxa over time
#'
#' Sums the relative abundance of the MAGs of each group in every sample and
#' draws its mean and standard deviation across replicates over time, with one
#' line per condition and one panel per group. It is the figure that answers
#' how a set of genomes (for example, those carrying hydrocarbon degradation
#' genes) responds to a disturbance.
#'
#' The same works without MAGs, on the profiles of the gene catalog or of the
#' reads: a functional profile from [read_feature_profile()] (KEGG orthologues
#' or hydrocarbon markers), or a taxonomic profile from
#' [read_taxonomy_profile()] or [read_singlem()] reduced to one rank. Then each
#' feature or taxon plays the role of a bin, and `groups` gathers them (with a
#' column `feature` or `taxon` instead of `bin`). Without `groups`, each one
#' gets its own panel, up to 16. A gene with several KOs of the same group
#' counts once per KO.
#'
#' @param abundance A data frame returned by [read_mag_abundance()],
#'   [read_feature_profile()], [read_taxonomy_profile()] or [read_singlem()].
#' @param samples A data frame with one row per sample: the column `sample`
#'   and the columns named by `time` and `condition`, for example from
#'   [parse_sample_names()].
#' @param groups Optional data frame with the columns `bin` and `group`, for
#'   example from [hydrocarbon_groups()]. Without it, a single panel shows the
#'   total abundance of all the bins in `abundance`. For a profile, the column
#'   `feature` or `taxon` replaces `bin`; [hydrocarbon_kos()] (`ko` and
#'   `pathway`) is also accepted.
#' @param time Column of `samples` with the sampling time. Numeric values are
#'   placed to scale; anything else is kept in order of appearance.
#' @param condition Optional column of `samples` used to draw one line per
#'   condition.
#' @param time_label Title of the x axis.
#' @param condition_labels Optional named character vector to show readable
#'   names in the legend instead of the codes used in the sample names, for
#'   example `c(C = "Low impact", O = "High impact")`.
#' @param condition_title Title of the legend.
#'
#' @return A `ggplot` object. The summarised values are in its `data` element.
#' @export
#'
#' @examples
#' dir <- system.file("extdata", "marine_honour", package = "magreport")
#' abundance <- read_mag_abundance(file.path(dir, "mag_abundance.tsv.gz"))
#' genes <- read_hydrocarbon_genes(file.path(dir, "hydrocarbon_per_mag.tsv"))
#' samples <- parse_sample_names(
#'   abundance$sample,
#'   into = c("time", "site", "condition", "replicate")
#' )
#' # Days after the oil arrived at each sampling time
#' days <- c(T0 = 2, T1 = 17, T2 = 33, T3 = 48, T4 = 63, T5 = 93, T6 = 124, T7 = 185)
#' samples$day <- days[samples$time]
#'
#' plot_abundance_timeline(
#'   abundance, samples,
#'   groups = hydrocarbon_groups(genes),
#'   time = "day", condition = "condition",
#'   time_label = "Days after oil arrival",
#'   condition_labels = c(C = "L (low)", O = "H (high)"),
#'   condition_title = "Oil impact"
#' )
plot_abundance_timeline <- function(abundance,
                                    samples,
                                    groups = NULL,
                                    time = "time",
                                    condition = NULL,
                                    time_label = time,
                                    condition_labels = NULL,
                                    condition_title = condition) {
  check_ggplot2()
  stopifnot(is.data.frame(abundance), is.data.frame(samples))

  # MAGs (bin), functional profile (feature) or taxonomic profile (taxon)
  unit <- intersect(c("bin", "feature", "taxon"), names(abundance))[1]
  if (is.na(unit) || !all(c("sample", "rel_abundance") %in% names(abundance))) {
    stop(
      "abundance must come from read_mag_abundance(), read_feature_profile(), ",
      "read_taxonomy_profile() or read_singlem()",
      call. = FALSE
    )
  }
  if (unit == "taxon") {
    if ("rank" %in% names(abundance) && length(unique(abundance$rank)) > 1) {
      stop(
        "abundance holds several ranks: keep one, for example ",
        "abundance[abundance$rank == \"genus\", ]",
        call. = FALSE
      )
    }
  }
  if (unit != "bin") {
    abundance$bin <- as.character(abundance[[unit]])
  }
  needed <- c("sample", time, condition)
  missing <- setdiff(needed, names(samples))
  if (length(missing) > 0) {
    stop("Column(s) not found in samples: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  absent <- setdiff(unique(abundance$sample), samples$sample)
  if (length(absent) > 0) {
    stop(
      length(absent), " sample(s) of abundance are missing from samples: ",
      paste(utils::head(absent, 3), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(groups)) {
    if (unit == "bin") {
      groups <- data.frame(bin = unique(abundance$bin), group = "All MAGs")
    } else {
      each <- unique(abundance$bin)
      if (length(each) > 16) {
        stop(
          length(each), " ", unit, "s in abundance: give groups, or keep 16 at most",
          call. = FALSE
        )
      }
      groups <- data.frame(bin = each, group = each, stringsAsFactors = FALSE)
    }
  } else if (!"bin" %in% names(groups)) {
    if (!"group" %in% names(groups) && "pathway" %in% names(groups)) {
      groups$group <- groups$pathway
    }
    id <- intersect(c(unit, "feature", "taxon", "ko", "marker"), names(groups))[1]
    if (!is.na(id)) {
      groups$bin <- as.character(groups[[id]])
    }
  }
  if (!all(c("bin", "group") %in% names(groups))) {
    stop("groups needs the columns ", unit, " and group", call. = FALSE)
  }

  # Sum per sample and group; samples where no MAG of a group was detected
  # count as zero, not as missing
  hits <- merge(abundance[c("bin", "sample", "rel_abundance")], groups[c("bin", "group")], by = "bin")
  if (nrow(hits) == 0) {
    stop("None of the ", unit, "s in groups is in abundance", call. = FALSE)
  }
  per_sample <- stats::aggregate(
    rel_abundance ~ sample + group, data = hits, FUN = sum
  )
  grid <- expand.grid(
    sample = unique(abundance$sample),
    group = unique(as.character(groups$group)),
    stringsAsFactors = FALSE
  )
  per_sample <- merge(grid, per_sample, by = c("sample", "group"), all.x = TRUE)
  per_sample$rel_abundance[is.na(per_sample$rel_abundance)] <- 0
  per_sample <- merge(per_sample, samples[needed], by = "sample")

  by <- c("group", time, condition)
  summary <- stats::aggregate(
    per_sample$rel_abundance,
    by = per_sample[by],
    FUN = function(v) c(mean = mean(v), sd = if (length(v) > 1) stats::sd(v) else 0, n = length(v))
  )
  summary <- cbind(summary[by], as.data.frame(summary$x))

  group_levels <- if (is.factor(groups$group)) levels(groups$group) else unique(groups$group)
  summary$group <- factor(summary$group, levels = group_levels[group_levels %in% summary$group])

  if (!is.numeric(summary[[time]])) {
    summary[[time]] <- factor(summary[[time]], levels = unique(samples[[time]]))
  }
  summary <- summary[order(summary$group, summary[[time]]), ]

  aes_line <- if (is.null(condition)) {
    ggplot2::aes(x = .data[[time]], y = .data$mean, group = 1)
  } else {
    ggplot2::aes(
      x = .data[[time]], y = .data$mean,
      colour = .data[[condition]], group = .data[[condition]]
    )
  }

  p <- ggplot2::ggplot(summary, aes_line) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(.data$mean - .data$sd, 0), ymax = .data$mean + .data$sd),
      width = 0, linewidth = 0.5, alpha = 0.7
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$group), scales = "free_y") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      x = time_label,
      y = "Relative abundance (%)",
      title = if (unit == "bin") "MAG abundance over time" else "Abundance over time",
      caption = paste(
        "Points: mean across replicates. Bars: standard deviation.",
        switch(unit,
          bin = "Relative abundance: share of the summed depth of all MAGs in the sample.",
          feature = "Relative abundance: share of the TPM of the gene catalog in the sample.",
          taxon = "Relative abundance: share of the sample at that rank."
        ),
        sep = "\n"
      )
    ) +
    theme_magreport()

  if (is.null(condition)) {
    p <- p + ggplot2::aes(colour = NULL) +
      ggplot2::scale_colour_manual(values = magreport_colours()[["olive"]], guide = "none")
  } else {
    if (is.null(condition_labels)) {
      p <- p + scale_colour_magreport(name = condition_title)
    } else {
      unknown <- setdiff(unique(as.character(summary[[condition]])), names(condition_labels))
      if (length(unknown) > 0) {
        stop("condition_labels has no label for: ", paste(unknown, collapse = ", "), call. = FALSE)
      }
      p <- p + scale_colour_magreport(name = condition_title, labels = condition_labels)
    }
  }

  p
}
