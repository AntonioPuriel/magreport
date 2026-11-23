#' Read the manifest of a mag-pipeline run
#'
#' Since version 2.1.0, mag-pipeline writes `manifest.json` in its results
#' directory: the pipeline version, the run, the databases with their dates,
#' every parameter and the list of summary tables the run produced. The format
#' is described in the
#' [pipeline documentation](https://github.com/AntonioPuriel/mag-pipeline/blob/main/docs/manifest.md).
#'
#' @param path `manifest.json`, or the directory that holds it.
#'
#' @return A list with the elements `schema_version`, `pipeline`, `run`,
#'   `databases`, `params` and `tables`. `tables` is a named list of character
#'   vectors with the full paths of the files of each table. The attribute
#'   `"dir"` holds the directory of the manifest.
#'
#'   An error is raised when the manifest has a major `schema_version` this
#'   version of magreport does not know, and a warning when the minor version
#'   is newer (the manifest may list tables magreport does not read yet).
#' @export
#'
#' @examples
#' dir <- tempfile()
#' dir.create(file.path(dir, "05_quality"), recursive = TRUE)
#' file.copy(
#'   system.file("extdata", "example_run", "mag_quality.tsv", package = "magreport"),
#'   file.path(dir, "05_quality")
#' )
#' writeLines(
#'   '{"schema_version": "1.0",
#'     "pipeline": {"name": "AntonioPuriel/mag-pipeline", "version": "2.1.0"},
#'     "tables": {"mag_quality": ["05_quality/mag_quality.tsv"]}}',
#'   file.path(dir, "manifest.json")
#' )
#' manifest <- read_manifest(dir)
#' manifest$pipeline$version
#' manifest$tables$mag_quality
read_manifest <- function(path) {
  if (dir.exists(path)) {
    path <- file.path(path, "manifest.json")
  }
  if (!file.exists(path)) {
    stop("Manifest not found: ", path, call. = FALSE)
  }

  manifest <- jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  check_schema_version(manifest$schema_version, path)

  dir <- dirname(normalizePath(path, winslash = "/"))
  tables <- manifest$tables
  if (is.null(tables) || length(tables) == 0) {
    tables <- list()
  }
  manifest$tables <- lapply(tables, function(files) {
    file.path(dir, unlist(files, use.names = FALSE))
  })

  for (name in c("pipeline", "run", "databases", "params")) {
    if (is.null(manifest[[name]])) {
      manifest[name] <- list(list())
    }
  }

  attr(manifest, "dir") <- dir
  manifest
}

# Highest schema_version of manifest.json this version of magreport reads
manifest_schema <- c(major = 1L, minor = 0L)

check_schema_version <- function(version, path) {
  if (is.null(version)) {
    stop("No schema_version in ", path, call. = FALSE)
  }
  parts <- suppressWarnings(as.integer(strsplit(as.character(version), ".", fixed = TRUE)[[1]]))
  if (length(parts) != 2 || anyNA(parts)) {
    stop("Invalid schema_version '", version, "' in ", path, call. = FALSE)
  }
  supported <- paste(manifest_schema, collapse = ".")
  if (parts[1] != manifest_schema[["major"]]) {
    stop(
      "The manifest has schema_version ", version, " and this version of ",
      "magreport reads ", supported, ". Install the version of magreport that ",
      "matches the pipeline (see docs/manifest.md in mag-pipeline).",
      call. = FALSE
    )
  }
  if (parts[2] > manifest_schema[["minor"]]) {
    warning(
      "The manifest has schema_version ", version, ", newer than the ",
      supported, " this version of magreport knows: some tables may be ",
      "left out. Update magreport.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# The tables of read_pipeline() and the manifest name of each one
manifest_table_names <- c(
  mags            = "mag_quality",
  bins            = "bin_summary",
  mapping         = "mapping_summary",
  depth           = "depth",
  assembly        = "assembly_stats",
  kos             = "ko_per_mag",
  annotation      = "annotation_stats",
  mimag           = "mimag_quality",
  hydrocarbon     = "hydrocarbon_per_mag",
  abundance       = "mag_abundance",
  catalog         = "gene_catalog_stats",
  profile_stats   = "profile_stats",
  contig_taxonomy = "taxonomy_per_sample",
  contig_kos      = "ko_per_sample",
  contig_canthyd  = "canthyd_per_sample",
  contig_hadeg    = "hadeg_per_sample"
)

# The files read_pipeline() uses, taken from the manifest. A table listed but
# missing on disk (for example, left behind when the results were copied) is
# reported and treated as absent.
manifest_files <- function(manifest) {
  tables <- manifest$tables
  missing <- character()

  existing <- function(name) {
    files <- tables[[name]]
    if (is.null(files)) {
      return(NULL)
    }
    gone <- files[!file.exists(files)]
    missing <<- c(missing, gone)
    files <- files[file.exists(files)]
    if (length(files) == 0) NULL else files
  }

  files <- lapply(manifest_table_names, function(name) {
    found <- existing(name)
    if (is.null(found)) NULL else found[1]
  })
  names(files) <- names(manifest_table_names)

  # SingleM: the profile of all samples, else those of each sample
  files["singlem"] <- list(existing("singlem")[1] %||% existing("singlem_per_sample"))
  # GUNC: one report per assembly and bin set, all of them
  files["gunc"] <- list(existing("gunc"))
  files["run_info"] <- list(NULL)

  if (length(missing) > 0) {
    warning(
      "The manifest lists ", length(missing), " file(s) that do not exist, ",
      "left out: ", paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }
  files
}
