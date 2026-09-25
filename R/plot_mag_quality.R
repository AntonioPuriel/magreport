#' Plot MAG completeness against contamination
#'
#' Draws the standard MAG quality plot: completeness on the x axis,
#' contamination on the y axis, with the MIMAG thresholds marked. Point size is
#' proportional to genome size when that column is available.
#'
#' @param mags A data frame returned by [read_mag_quality()].
#' @param facet_by Optional column used to split the plot into panels, for
#'   example `"binner"`.
#' @param label_high If `TRUE`, writes the bin name next to high-quality MAGs.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "mag_quality.tsv",
#'                     package = "magreport")
#' mags <- read_mag_quality(path)
#' plot_mag_quality(mags)
#' plot_mag_quality(mags, facet_by = "binner")
plot_mag_quality <- function(mags, facet_by = NULL, label_high = FALSE) {
  check_ggplot2()
  stopifnot(is.data.frame(mags))

  required <- c("completeness", "contamination")
  missing <- setdiff(required, names(mags))
  if (length(missing) > 0) {
    stop("Column(s) not found in mags: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  mags$quality <- factor(mags$quality, levels = quality_levels())

  p <- ggplot2::ggplot(
    mags,
    ggplot2::aes(x = .data$completeness, y = .data$contamination)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = 90, xmax = 100, ymin = 0, ymax = 5,
      fill = quality_colours()[["high"]], alpha = 0.12
    ) +
    ggplot2::geom_vline(
      xintercept = c(50, 90), linetype = "dashed", colour = magreport_colours()[["sand_mid"]]
    ) +
    ggplot2::geom_hline(
      yintercept = c(5, 10), linetype = "dashed", colour = magreport_colours()[["sand_mid"]]
    )

  if ("total_length" %in% names(mags)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = .data$quality, size = .data$total_length / 1e6),
      alpha = 0.85
    ) +
      ggplot2::scale_size_continuous(name = "Size (Mb)", range = c(1.5, 6))
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$quality), size = 3, alpha = 0.85)
  }

  if (isTRUE(label_high) && "bin" %in% names(mags)) {
    high <- mags[!is.na(mags$quality) & mags$quality == "high", , drop = FALSE]
    if (nrow(high) > 0) {
      p <- p + ggplot2::geom_text(
        data = high,
        ggplot2::aes(label = .data$bin),
        hjust = 1.1, vjust = -0.6, size = 3, show.legend = FALSE,
        colour = magreport_colours()[["olive_dark"]]
      )
    }
  }

  if (!is.null(facet_by)) {
    if (!facet_by %in% names(mags)) {
      stop("Column not found in mags: ", facet_by, call. = FALSE)
    }
    p <- p + ggplot2::facet_wrap(facet_by)
  }

  p +
    ggplot2::scale_colour_manual(
      name = "Quality", values = quality_colours(), drop = FALSE
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 100)) +
    ggplot2::labs(
      x = "Completeness (%)",
      y = "Contamination (%)",
      title = "MAG quality",
      caption = "Dashed lines: MIMAG thresholds"
    ) +
    theme_magreport()
}
