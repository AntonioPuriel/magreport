#' Plot the abundance of groups of MAGs over time
#'
#' Sums the relative abundance of the MAGs of each group in every sample and
#' draws its mean and standard deviation across replicates over time, with one
#' line per condition and one panel per group. It is the figure that answers
#' how a set of genomes (for example, those carrying hydrocarbon degradation
#' genes) responds to a disturbance.
#'
#' @param abundance A data frame returned by [read_mag_abundance()].
#' @param samples A data frame with one row per sample: the column `sample`
#'   and the columns named by `time` and `condition`, for example from
#'   [parse_sample_names()].
#' @param groups Optional data frame with the columns `bin` and `group`, for
#'   example from [hydrocarbon_groups()]. Without it, a single panel shows the
#'   total abundance of all the bins in `abundance`.
#' @param time Column of `samples` with the sampling time. Numeric values are
#'   placed to scale; anything else is kept in order of appearance.
#' @param condition Optional column of `samples` used to draw one line per
#'   condition.
#' @param time_label Title of the x axis.
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
#'   time_label = "Days after oil arrival"
#' )
plot_abundance_timeline <- function(abundance,
                                    samples,
                                    groups = NULL,
                                    time = "time",
                                    condition = NULL,
                                    time_label = time) {
  check_ggplot2()
  stopifnot(is.data.frame(abundance), is.data.frame(samples))

  if (!all(c("bin", "sample", "rel_abundance") %in% names(abundance))) {
    stop("abundance must come from read_mag_abundance()", call. = FALSE)
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
    groups <- data.frame(bin = unique(abundance$bin), group = "All MAGs")
  }
  if (!all(c("bin", "group") %in% names(groups))) {
    stop("groups needs the columns bin and group", call. = FALSE)
  }

  # Sum per sample and group; samples where no MAG of a group was detected
  # count as zero, not as missing
  hits <- merge(abundance[c("bin", "sample", "rel_abundance")], groups[c("bin", "group")], by = "bin")
  if (nrow(hits) == 0) {
    stop("None of the bins in groups is in abundance", call. = FALSE)
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
      title = "MAG abundance over time",
      caption = paste(
        "Points: mean across replicates. Bars: standard deviation.",
        "Relative abundance: share of the summed depth of all MAGs in the sample.",
        sep = "\n"
      )
    ) +
    theme_magreport()

  if (is.null(condition)) {
    p <- p + ggplot2::aes(colour = NULL) +
      ggplot2::scale_colour_manual(values = magreport_colours()[["olive"]], guide = "none")
  } else {
    p <- p + scale_colour_magreport(name = condition)
  }

  p
}
