#' Read the abundance of each MAG in each sample
#'
#' Reads `04_binning/mag_abundance.tsv`, the mean read depth of every bin in
#' every sample, and adds its relative abundance: the share of the summed depth
#' of all the bins kept in that sample.
#'
#' A bin that DAS Tool kept unchanged is listed twice in the pipeline table
#' (under its binner and under `dastool`) with the same depths: it is kept
#' once.
#'
#' @param path Path to `mag_abundance.tsv`.
#' @param bins Optional character vector of bins to keep, for example the
#'   dereplicated representatives (`read_mimag_quality(...)$bin`). The relative
#'   abundance is computed over these bins only.
#'
#' @return A data frame with one row per bin and sample, with the columns
#'   `bin`, `sample`, `mean_depth` and `rel_abundance` (in percent), plus
#'   `assembly`, `binner` and `bin_length` when present.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "marine_honour", "mag_abundance.tsv.gz",
#'                     package = "magreport")
#' abundance <- read_mag_abundance(path)
#' head(abundance)
read_mag_abundance <- function(path, bins = NULL) {
  x <- read_pipeline_table(path, required = c("bin", "sample", "mean_depth"))
  x$mean_depth <- as.numeric(x$mean_depth)

  if (!is.null(bins)) {
    missing <- setdiff(bins, x$bin)
    if (length(missing) > 0) {
      warning(
        length(missing), " bin(s) not found in ", basename(path), ": ",
        paste(utils::head(missing, 3), collapse = ", "),
        if (length(missing) > 3) ", ..." else "",
        call. = FALSE
      )
    }
    x <- x[x$bin %in% bins, , drop = FALSE]
  }

  x <- x[!duplicated(x[c("bin", "sample")]), , drop = FALSE]

  totals <- tapply(x$mean_depth, x$sample, sum)
  x$rel_abundance <- ifelse(
    totals[x$sample] > 0,
    100 * x$mean_depth / totals[x$sample],
    0
  )

  rownames(x) <- NULL
  x
}

#' Split sample names into metadata columns
#'
#' Many sample names encode the design, for example `T3-BB-O-2` for time T3,
#' site BB, condition O and replicate 2. This splits them into columns so they
#' can be used as the `samples` table of [plot_abundance_timeline()].
#'
#' @param samples Character vector of sample names.
#' @param into Names of the columns, one per part of the name.
#' @param sep Regular expression separating the parts.
#'
#' @return A data frame with the column `sample` and one column per element of
#'   `into`.
#' @export
#'
#' @examples
#' parse_sample_names(
#'   c("T0-BB-C-1", "T0-BB-O-1"),
#'   into = c("time", "site", "condition", "replicate")
#' )
parse_sample_names <- function(samples, into, sep = "-") {
  samples <- unique(as.character(samples))
  parts <- strsplit(samples, sep)

  wrong <- lengths(parts) != length(into)
  if (any(wrong)) {
    stop(
      "Sample names do not have ", length(into), " parts separated by '", sep,
      "': ", paste(utils::head(samples[wrong], 3), collapse = ", "),
      call. = FALSE
    )
  }

  out <- as.data.frame(do.call(rbind, parts), stringsAsFactors = FALSE)
  names(out) <- into
  cbind(data.frame(sample = samples, stringsAsFactors = FALSE), out)
}
