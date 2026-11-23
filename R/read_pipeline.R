#' Find and read every table of a mag-pipeline run
#'
#' Reads the output tables of mag-pipeline from a directory. When the
#' directory holds `manifest.json` (mag-pipeline 2.1.0 or later), the tables
#' it lists are read ([read_manifest()]). Otherwise the tables are looked for by
#' name, whatever the folder layout is: the `results/` directory of a run, or a
#' folder where the tables were copied together.
#'
#' @param results_dir Directory of the run: with `manifest.json` at its top, or
#'   searched recursively.
#'
#' @return A list with the elements `mags`, `mimag`, `bins`, `mapping`,
#'   `depth`, `assembly`, `kos`, `annotation`, `hydrocarbon` and `abundance`,
#'   each one a data frame or `NULL` when the table is not found. When the MIMAG
#'   table exists, `abundance` is restricted to its genomes (the dereplicated
#'   representatives).
#'
#'   The profiles of the community without MAGs (mag-pipeline 2.0.0 or later)
#'   are in `catalog` ([read_catalog_stats()]), `profile_stats`
#'   ([read_profile_stats()]), `contig_taxonomy` ([read_taxonomy_profile()]),
#'   `contig_kos`, `contig_canthyd`, `contig_hadeg` ([read_feature_profile()])
#'   and `singlem` ([read_singlem()]). For SingleM, the table of all samples is
#'   preferred, then the profile of each sample, then the tables per rank.
#'   The attribute `"files"` holds the paths that were used, and the attribute
#'   `"manifest"` the manifest of the run (`NULL` without one).
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

  manifest <- NULL
  if (file.exists(file.path(results_dir, "manifest.json"))) {
    manifest <- read_manifest(results_dir)
    files <- manifest_files(manifest)
  } else {
    files <- find_pipeline_files(results_dir)
  }

  out <- read_pipeline_tables(files)
  attr(out, "files") <- files
  attr(out, "manifest") <- manifest
  out
}

# The tables of a run without manifest.json (mag-pipeline 2.0.0 or earlier, or
# tables copied together), found by name
find_pipeline_files <- function(results_dir) {
  find_one <- function(pattern) {
    hits <- list.files(results_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) NULL else hits[1]
  }

  list(
    mags     = find_one("^mag_quality\\.tsv$"),
    bins     = find_one("^bin_summary\\.tsv$"),
    mapping  = find_one("^mapping_summary\\.tsv$"),
    depth    = find_one("\\.depth\\.txt$"),
    assembly = find_one("\\.contigs\\.stats\\.tsv$"),
    kos        = find_one("^ko_per_mag\\.tsv$"),
    annotation = find_one("^annotation_stats\\.tsv$"),
    mimag       = find_one("^mimag_quality\\.tsv$"),
    hydrocarbon = find_one("^hydrocarbon_per_mag\\.tsv$"),
    abundance   = find_one("^mag_abundance\\.tsv(\\.gz)?$"),
    catalog         = find_one("^gene_catalog\\.stats\\.tsv$"),
    profile_stats   = find_one("^profile_stats\\.tsv$"),
    contig_taxonomy = find_one("^taxonomy_per_sample\\.tsv(\\.gz)?$"),
    contig_kos      = find_one("^ko_per_sample\\.tsv(\\.gz)?$"),
    contig_canthyd  = find_one("^canthyd_per_sample\\.tsv(\\.gz)?$"),
    contig_hadeg    = find_one("^hadeg_per_sample\\.tsv(\\.gz)?$"),
    singlem         = find_singlem(results_dir)
  )
}

read_pipeline_tables <- function(files) {
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

  out["catalog"] <- list(read_if(files$catalog, read_catalog_stats))
  out["profile_stats"] <- list(read_if(files$profile_stats, read_profile_stats))
  out["contig_taxonomy"] <- list(read_if(files$contig_taxonomy, read_taxonomy_profile))
  for (name in c("contig_kos", "contig_canthyd", "contig_hadeg")) {
    out[name] <- list(read_if(files[[name]], function(p) {
      read_feature_profile(p, total_tpm = out$profile_stats)
    }))
  }
  out["singlem"] <- list(read_if(files$singlem, read_singlem_tables))
  out
}

# The SingleM tables of a run: the profile of all samples, else those of each
# sample, else the taxon-by-sample tables of every rank
find_singlem <- function(results_dir) {
  find_all <- function(pattern) {
    sort(list.files(results_dir, pattern = pattern, recursive = TRUE, full.names = TRUE))
  }
  extras <- find_all("^singlem_profiles_with_extras\\.tsv(\\.gz)?$")
  if (length(extras) > 0) {
    return(extras[1])
  }
  per_sample <- find_all("\\.singlem\\.tsv(\\.gz)?$")
  if (length(per_sample) > 0) {
    return(per_sample)
  }
  ranks <- paste(taxonomic_ranks(), collapse = "|")
  by_rank <- find_all(paste0("^singlem_(", ranks, ")\\.tsv(\\.gz)?$"))
  if (length(by_rank) > 0) by_rank else NULL
}

read_singlem_tables <- function(paths) {
  if (length(paths) > 1 && all(grepl("^singlem_[a-z]+\\.tsv", basename(paths)))) {
    out <- do.call(rbind, lapply(paths, read_singlem))
    attr(out, "source") <- "SingleM"
    out
  } else {
    read_singlem(paths)
  }
}
