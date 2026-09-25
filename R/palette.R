#' Colours of the magreport palette
#'
#' The key colours used across magreport: olive greens and desert tones. Every
#' plot of the package takes its colours from here, so they can also be used to
#' match other figures or slides to the report.
#'
#' @return A named character vector of hex colours.
#' @export
#' @seealso [magreport_palette()] for categorical data, [theme_magreport()].
#' @examples
#' magreport_colours()
#' magreport_colours()[["olive"]]
magreport_colours <- function() {
  c(
    olive       = "#5E6B2F",
    olive_dark  = "#353D1F",
    sand        = "#F3EBDD",
    sand_mid    = "#D8C8A8",
    taupe       = "#9C8B6E",
    terracotta  = "#A8552F",
    ochre       = "#D4A24C"
  )
}

# Categorical colours, ordered so that neighbouring levels contrast
magreport_discrete <- c(
  "#5E6B2F", # olive
  "#A8552F", # terracotta
  "#D4A24C", # ochre
  "#5F7F7A", # sage blue
  "#353D1F", # dark olive
  "#C98B6B", # clay
  "#A89F5B", # khaki
  "#7A3B1E", # rust
  "#8FA37A", # sage
  "#8C6A74", # mauve
  "#C9A57A", # tan
  "#9C8B6E"  # taupe
)

# Sequential ramp, from empty to full
magreport_sequential <- c("#F3EBDD", "#C9B98A", "#5E6B2F", "#2B3119")

#' Categorical colours for any number of groups
#'
#' Returns up to 12 distinct colours of the magreport palette. For more groups
#' the colours are interpolated, which makes neighbouring ones harder to tell
#' apart: consider grouping the smallest categories first.
#'
#' @param n Number of colours.
#'
#' @return A character vector of `n` hex colours.
#' @export
#' @examples
#' magreport_palette(4)
magreport_palette <- function(n) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 0)
  n <- as.integer(n)
  if (n <= length(magreport_discrete)) {
    return(magreport_discrete[seq_len(n)])
  }
  grDevices::colorRampPalette(magreport_discrete)(n)
}

# discrete_scale() dropped its scale_name argument in ggplot2 3.5.0
magreport_discrete_scale <- function(aesthetics, ...) {
  check_ggplot2()
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::discrete_scale(aesthetics, palette = magreport_palette, ...)
  } else {
    ggplot2::discrete_scale(aesthetics, "magreport", palette = magreport_palette, ...)
  }
}

#' Discrete ggplot2 scales with the magreport palette
#'
#' @param ... Passed to [ggplot2::discrete_scale()], for example `name`.
#'
#' @return A ggplot2 scale.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, colour = factor(cyl))) +
#'     ggplot2::geom_point(size = 3) +
#'     scale_colour_magreport(name = "Cylinders")
#' }
scale_colour_magreport <- function(...) {
  magreport_discrete_scale("colour", ...)
}

#' @rdname scale_colour_magreport
#' @export
scale_color_magreport <- scale_colour_magreport

#' @rdname scale_colour_magreport
#' @export
scale_fill_magreport <- function(...) {
  magreport_discrete_scale("fill", ...)
}

#' Continuous fill scale with the magreport palette
#'
#' A sequential scale from sand to dark olive, used for heatmaps.
#'
#' @param ... Passed to [ggplot2::scale_fill_gradientn()], for example `name`.
#'
#' @return A ggplot2 scale.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(faithfuld, ggplot2::aes(waiting, eruptions, fill = density)) +
#'     ggplot2::geom_raster() +
#'     scale_fill_magreport_c()
#' }
scale_fill_magreport_c <- function(...) {
  check_ggplot2()
  ggplot2::scale_fill_gradientn(colours = magreport_sequential, ...)
}

#' The ggplot2 theme of magreport
#'
#' A light theme based on [ggplot2::theme_bw()] with the olive and sand tones
#' of the package. Every magreport plot uses it, and it can be added to any
#' other ggplot.
#'
#' @param base_size Base font size, in points.
#'
#' @return A ggplot2 theme.
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'     ggplot2::geom_point(colour = magreport_colours()[["olive"]]) +
#'     theme_magreport()
#' }
theme_magreport <- function(base_size = 11) {
  check_ggplot2()
  col <- magreport_colours()
  line <- "#DDD3C1"

  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = col[["olive_dark"]]),
      plot.title = ggplot2::element_text(face = "bold", colour = col[["olive_dark"]]),
      plot.caption = ggplot2::element_text(colour = col[["taupe"]]),
      axis.text = ggplot2::element_text(colour = "#4A4F3A"),
      axis.ticks = ggplot2::element_line(colour = line),
      panel.border = ggplot2::element_rect(colour = line, fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "#EFE9DE"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = col[["sand"]], colour = line),
      strip.text = ggplot2::element_text(colour = col[["olive_dark"]], face = "bold"),
      legend.key = ggplot2::element_blank()
    )
}
