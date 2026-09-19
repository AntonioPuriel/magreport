#' @keywords internal
read_pipeline_table <- function(path, required = character()) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  x <- utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    # GTDB-Tk writes N/A, CheckM2 and the pipeline write NA
    na.strings = c("NA", "N/A", "na", "")
  )

  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "Unexpected format in ", basename(path), ": missing column(s) ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x
}

#' Quality categories used by magreport
#'
#' The levels follow the MIMAG standard, from the best to the worst category.
#'
#' @return A character vector with the quality levels.
#' @export
#' @examples
#' quality_levels()
quality_levels <- function() {
  c("high", "medium", "low")
}

#' Colours used for the quality categories
#'
#' @return A named character vector of colours, one per level of
#'   [quality_levels()].
#' @export
#' @examples
#' quality_colours()
quality_colours <- function() {
  c(high = "#1B9E77", medium = "#D95F02", low = "#7570B3")
}

#' @keywords internal
check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is needed for the plotting functions. ",
      "Install it with install.packages(\'ggplot2\').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
