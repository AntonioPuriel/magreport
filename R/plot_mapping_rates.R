#' Plot the alignment rate of each sample
#'
#' A low alignment rate means that a large part of a sample is not represented
#' in the assembly, which is worth checking before interpreting coverage. When
#' samples were mapped against several assemblies (mag-pipeline with
#' `--assembly_mode group` or `per_sample`), each assembly gets its own panel.
#'
#' @param mapping A data frame returned by [read_mapping_summary()].
#' @param warn_below Alignment rate under which bars are highlighted, in
#'   percent. Use `NULL` to disable the highlight.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "mapping_summary.tsv",
#'                     package = "magreport")
#' mapping <- read_mapping_summary(path)
#' plot_mapping_rates(mapping)
plot_mapping_rates <- function(mapping, warn_below = 70) {
  check_ggplot2()
  stopifnot(is.data.frame(mapping))

  if (!"overall_alignment_rate_pct" %in% names(mapping)) {
    stop("Column not found in mapping: overall_alignment_rate_pct", call. = FALSE)
  }

  # Order samples by their mean alignment rate; a sample can appear once per assembly
  mean_rate <- tapply(mapping$overall_alignment_rate_pct, mapping$sample, mean, na.rm = TRUE)
  mapping$sample <- factor(mapping$sample, levels = names(sort(mean_rate)))
  several_assemblies <- "assembly" %in% names(mapping) && length(unique(mapping$assembly)) > 1
  mapping$flag <- if (is.null(warn_below)) {
    "ok"
  } else {
    ifelse(mapping$overall_alignment_rate_pct < warn_below, "low", "ok")
  }

  p <- ggplot2::ggplot(
    mapping,
    ggplot2::aes(x = .data$sample, y = .data$overall_alignment_rate_pct, fill = .data$flag)
  ) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", .data$overall_alignment_rate_pct)),
      hjust = -0.15, size = 3
    ) +
    ggplot2::scale_fill_manual(values = c(ok = "#2C7FB8", low = "#D95F02")) +
    ggplot2::scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Reads aligned to the assembly (%)",
      title = "Read mapping per sample"
    ) +
    ggplot2::theme_bw()

  if (several_assemblies) {
    p <- p + ggplot2::facet_wrap(~assembly)
  }

  p
}
