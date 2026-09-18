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
    na.strings = c("NA", "")
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
