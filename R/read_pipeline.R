#' Find and read every table of a mag-pipeline run
#'
#' Looks for the output tables of mag-pipeline inside a directory, whatever the
#' folder layout is: the `results/` directory of a run, or a folder where the
#' tables were copied together.
#'
#' @param results_dir Directory to search, recursively.
#'
#' @return A list with the elements `mags`, `mimag`, `bins`, `mapping`,
#'   `depth`, `assembly`, `kos`, `annotation`, `hydrocarbon` and `abundance`,
#'   each one a data frame or `NULL` when the table is not found. When the MIMAG
#'   table exists, `abundance` is restricted to its genomes (the dereplicated
#'   representatives).
#'   The attribute `"files"` holds the paths that were used.
#' @export
#'
#' @examples
#' results_dir <- system.file("extdata", "example_run", package = "magreport")
#' run <- read_pipeline(results_dir)
#' names(run)
#' nrow(run$mags)
read_pipeline <- function(results_dir) {
  if (!dir.exists(results_dir)) {
    stop("Directory not found: ", results_dir, call. = FALSE)
  }

  find_one <- function(pattern) {
    hits <- list.files(results_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) NULL else hits[1]
  }

  files <- list(
    mags     = find_one("^mag_quality\\.tsv$"),
    bins     = find_one("^bin_summary\\.tsv$"),
    mapping  = find_one("^mapping_summary\\.tsv$"),
    depth    = find_one("\\.depth\\.txt$"),
    assembly = find_one("\\.contigs\\.stats\\.tsv$"),
    kos        = find_one("^ko_per_mag\\.tsv$"),
    annotation = find_one("^annotation_stats\\.tsv$"),
    mimag       = find_one("^mimag_quality\\.tsv$"),
    hydrocarbon = find_one("^hydrocarbon_per_mag\\.tsv$"),
    abundance   = find_one("^mag_abundance\\.tsv(\\.gz)?$")
  )

  read_if <- function(path, reader) {
    if (is.null(path)) NULL else reader(path)
  }

  out <- list(
    mags     = read_if(files$mags, read_mag_quality),
    bins     = read_if(files$bins, read_bin_summary),
    mapping  = read_if(files$mapping, read_mapping_summary),
    depth    = read_if(files$depth, read_depth),
    assembly = read_if(files$assembly, function(p) read_pipeline_table(p)),
    kos        = read_if(files$kos, read_ko_table),
    annotation = read_if(files$annotation, read_annotation_stats),
    mimag       = read_if(files$mimag, read_mimag_quality),
    hydrocarbon = read_if(files$hydrocarbon, read_hydrocarbon_genes)
  )

  # list() keeps the element when the table is missing (NULL)
  out["abundance"] <- list(read_if(files$abundance, function(p) {
    read_mag_abundance(p, bins = if (is.null(out$mimag)) NULL else out$mimag$bin)
  }))

  attr(out, "files") <- files
  out
}
