#' Read the contig depth table
#'
#' Reads the table written by `jgi_summarize_bam_contig_depths`
#' (`03_mapping/<assembly>.depth.txt`), which holds the coverage of every
#' contig in every sample.
#'
#' @param path Path to the depth file.
#' @param long If `TRUE` (the default), returns one row per contig and sample.
#'   If `FALSE`, returns the table as written by the tool.
#'
#' @return A data frame. In long format it has the columns `contig`,
#'   `contig_length`, `sample`, `depth` and `depth_var`.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_run", "coassembly.depth.txt",
#'                     package = "magreport")
#' depth <- read_depth(path)
#' head(depth)
read_depth <- function(path, long = TRUE) {
  x <- read_pipeline_table(
    path,
    required = c("contigName", "contigLen", "totalAvgDepth")
  )

  if (!long) {
    return(x)
  }

  # Everything after the fixed columns is per-sample: <name> and <name>-var
  # (read.delim turns the dash into a dot). Some versions of the tool also
  # write the columns without the .bam suffix.
  sample_cols <- setdiff(names(x), c("contigName", "contigLen", "totalAvgDepth"))
  var_cols <- grep("[.-]var$", sample_cols, value = TRUE)
  depth_cols <- setdiff(sample_cols, var_cols)

  if (length(depth_cols) == 0) {
    stop("No per-sample depth columns found in ", basename(path), call. = FALSE)
  }

  pieces <- lapply(depth_cols, function(column) {
    sample <- sub("\\.bam$", "", column)
    sample <- sub("\\.sorted$", "", sample)

    var_col <- var_cols[sub("[.-]var$", "", var_cols) == column]

    data.frame(
      contig = x$contigName,
      contig_length = as.numeric(x$contigLen),
      sample = sample,
      depth = as.numeric(x[[column]]),
      depth_var = if (length(var_col) > 0) as.numeric(x[[var_col[1]]]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}
