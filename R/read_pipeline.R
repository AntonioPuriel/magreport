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
#' @param depth If `FALSE`, the contig depth tables are not read. They have one
#'   row per contig and sample, several GB in large runs, and [mag_report()]
#'   does not use them.
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
#'   and `singlem` ([read_singlem()]). The GUNC reports of every bin set
#'   (`--gunc_db`) are in `gunc` ([read_gunc()]). For SingleM, the table of all samples is
#'   preferred, then the profile of each sample, then the tables per rank.
#'   The attribute `"files"` holds the paths that were used, and the attribute
#'   `"manifest"` the manifest of the run (`NULL` without one). The attribute
#'   `"run_info"` describes the run (pipeline version, databases, parameters):
#'   the manifest, or else `run_info.txt` read with [read_run_info()].
#'
#'   **Several runs.** With several directories, or a named one, every table
#'   of every run is stacked with a first column `study` (the names of
#'   `results_dir`; without names, the directory names). Columns that only
#'   some runs have are filled with `NA`. Bin and sample names may repeat
#'   between studies, so join tables by `study` and `bin`, as [plot_gunc()]
#'   and [mag_recovery()] do. The attributes `"files"`, `"manifest"` and
#'   `"run_info"` become lists with one element per study, and the attribute
#'   `"studies"` holds the study names.
#' @export
#'
#' @examples
#' results_dir <- system.file("extdata", "example_run", package = "magreport")
#' run <- read_pipeline(results_dir)
#' names(run)
#' nrow(run$mags)
#'
#' # Several runs, for example the studies of a catalog
#' runs <- read_pipeline(c(
#'   gut = system.file("extdata", "example_run", package = "magreport"),
#'   sediment = system.file("extdata", "marine_honour", package = "magreport")
#' ))
#' table(runs$mags$study)
#' table(runs$mimag$study)
read_pipeline <- function(results_dir, depth = TRUE) {
  if (length(results_dir) == 0) {
    stop("No results directory given", call. = FALSE)
  }
  if (length(results_dir) == 1 && is.null(names(results_dir))) {
    return(read_pipeline_one(results_dir, depth = depth))
  }

  studies <- names(results_dir)
  if (is.null(studies)) {
    studies <- basename(normalizePath(results_dir, mustWork = FALSE))
  }
  if (anyNA(studies) || any(studies == "")) {
    stop("Every results directory needs a study name", call. = FALSE)
  }
  if (anyDuplicated(studies)) {
    stop("Study names must be unique: ", paste(unique(studies[duplicated(studies)]), collapse = ", "),
         "; name them, as in c(study1 = dir1, study2 = dir2)", call. = FALSE)
  }

  runs <- lapply(results_dir, read_pipeline_one, depth = depth)
  names(runs) <- studies

  out <- lapply(stats::setNames(names(runs[[1]]), names(runs[[1]])), function(table) {
    bind_studies(lapply(runs, `[[`, table))
  })
  for (what in c("files", "manifest", "run_info")) {
    attr(out, what) <- lapply(runs, attr, what)
  }
  attr(out, "studies") <- studies
  out
}

# Stack the same table of several studies, with a first column `study`.
# Columns missing in some studies are filled with NA; factors stay factors
# when their levels agree.
bind_studies <- function(tables) {
  tables <- Filter(Negate(is.null), tables)
  if (length(tables) == 0) {
    return(NULL)
  }
  columns <- unique(unlist(lapply(tables, names)))
  filled <- lapply(names(tables), function(study) {
    x <- tables[[study]]
    for (column in setdiff(columns, names(x))) {
      x[[column]] <- rep(NA, nrow(x))
    }
    x <- x[, columns, drop = FALSE]
    cbind(data.frame(study = rep(study, nrow(x)), stringsAsFactors = FALSE), x)
  })
  out <- do.call(rbind, filled)
  rownames(out) <- NULL
  source <- unique(unlist(lapply(tables, attr, "source")))
  if (length(source) > 0) {
    attr(out, "source") <- paste(source, collapse = " + ")
  }
  out
}

# Join key of a table: the bin, within its study when there are several
bin_key <- function(x) {
  if ("study" %in% names(x)) paste(x$study, x$bin, sep = "\r") else x$bin
}

read_pipeline_one <- function(results_dir, depth = TRUE) {
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

  if (!isTRUE(depth)) {
    files["depth"] <- list(NULL)
  }
  out <- read_pipeline_tables(files)
  attr(out, "files") <- files
  attr(out, "manifest") <- manifest
  attr(out, "run_info") <- if (!is.null(manifest)) {
    manifest
  } else if (!is.null(files$run_info)) {
    read_run_info(files$run_info)
  }
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
    depth    = find_all(results_dir, "\\.depth\\.txt(\\.gz)?$"),
    assembly = find_all(results_dir, "\\.contigs\\.stats\\.tsv(\\.gz)?$"),
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
    singlem         = find_singlem(results_dir),
    gunc            = find_all(results_dir, "\\.gunc\\.tsv(\\.gz)?$"),
    run_info        = find_one("^run_info\\.txt$")
  )
}

find_all <- function(results_dir, pattern) {
  hits <- sort(list.files(results_dir, pattern = pattern, recursive = TRUE, full.names = TRUE))
  if (length(hits) == 0) NULL else hits
}

read_pipeline_tables <- function(files) {
  read_if <- function(path, reader) {
    if (is.null(path)) NULL else reader(path)
  }

  out <- list(
    mags     = read_if(files$mags, read_mag_quality),
    bins     = read_if(files$bins, read_bin_summary),
    mapping  = read_if(files$mapping, read_mapping_summary),
    depth    = read_if(files$depth, read_depths),
    assembly = read_if(files$assembly, function(paths) {
      do.call(rbind, lapply(paths, read_pipeline_table))
    }),
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
  out["gunc"] <- list(read_if(files$gunc, read_gunc))
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

#' Summarise several runs side by side
#'
#' One row per study of runs read together with [read_pipeline()]: samples,
#' bins, MAGs of each quality, bins that fail GUNC and pipeline version. The
#' MAGs are those of the MIMAG table (the dereplicated genomes) when a study
#' has one, else those of the CheckM2 table.
#'
#' @param run A list returned by [read_pipeline()] with several directories.
#'
#' @return A data frame with one row per study and the columns `study`,
#'   `pipeline`, `samples`, `bins`, `high`, `medium`, `low` and `gunc_fail`
#'   (`NA` when the table behind a column is missing).
#' @export
#'
#' @examples
#' runs <- read_pipeline(c(
#'   gut = system.file("extdata", "example_run", package = "magreport"),
#'   sediment = system.file("extdata", "marine_honour", package = "magreport")
#' ))
#' summarise_studies(runs)
summarise_studies <- function(run) {
  studies <- attr(run, "studies")
  if (is.null(studies)) {
    stop("run holds a single run: read several with read_pipeline(c(a = dir1, b = dir2))",
         call. = FALSE)
  }

  in_study <- function(x, study) {
    if (is.null(x)) NULL else x[x$study == study, , drop = FALSE]
  }

  rows <- lapply(studies, function(study) {
    samples <- unique(c(in_study(run$mapping, study)$sample,
                        in_study(run$singlem, study)$sample,
                        in_study(run$profile_stats, study)$sample))
    mimag <- in_study(run$mimag, study)
    mags <- if (!is.null(mimag) && nrow(mimag) > 0) mimag else in_study(run$mags, study)
    count <- function(level) {
      if (is.null(mags) || nrow(mags) == 0) NA_integer_ else sum(mags$quality == level, na.rm = TRUE)
    }
    gunc <- in_study(run$gunc, study)
    info <- attr(run, "run_info")[[study]]
    data.frame(
      study = study,
      pipeline = if (is.null(info$pipeline$version)) NA_character_ else info$pipeline$version,
      samples = if (length(samples) == 0) NA_integer_ else length(samples),
      bins = if (is.null(run$bins) || !study %in% run$bins$study) NA_integer_ else
        nrow(in_study(run$bins, study)),
      high = count("high"),
      medium = count("medium"),
      low = count("low"),
      gunc_fail = if (is.null(gunc) || nrow(gunc) == 0) NA_integer_ else sum(!gunc$pass, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# The depth tables of every assembly, stacked, with a first column `assembly`
# taken from the file name (<assembly>.depth.txt)
read_depths <- function(paths) {
  tables <- lapply(paths, function(path) {
    x <- read_depth(path)
    assembly <- sub("\\.depth\\.txt(\\.gz)?$", "", basename(path))
    cbind(data.frame(assembly = rep(assembly, nrow(x)), stringsAsFactors = FALSE), x)
  })
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}
