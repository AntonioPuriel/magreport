#' Write every table of a run, ready for other tools
#'
#' Writes the results of a mag-pipeline run as plain tables, so they can be
#' analysed with any other tool (R, Python, Excel, vegan, phyloseq...) without
#' reading the pipeline output again:
#'
#' * **Summaries** (one row per item): `samples.tsv`, `assemblies.tsv`,
#'   `bins.tsv` (every bin, with CheckM2, GTDB-Tk, GUNC and annotation) and
#'   `mags.tsv` (the dereplicated MAGs, with their abundance and functions).
#' * **Matrices** (`matrices/`, one row per feature and one column per sample):
#'   the abundance of the MAGs, the CAT and SingleM profiles of every rank and
#'   the KEGG, CANT-HYD and HADEG profiles of the gene catalog, each one twice:
#'   as the measure of the pipeline (`depth`, `tpm`, `coverage`) and as a
#'   percentage of the sample (`relative`). Also the KOs and hydrocarbon
#'   markers of each MAG.
#' * **Long tables** (`long/`): the tables of [read_pipeline()] as they are,
#'   one row per feature and sample, for ggplot2 or a database.
#' * **Methods** (`methods/`): pipeline version, databases, tool versions and
#'   parameters ([describe_run()]).
#' * **Excel** (`magreport_tables.xlsx`, needs the package writexl): the
#'   summaries, the methods and every matrix, one sheet each.
#' * **phyloseq** (`phyloseq/`, needs the package phyloseq): `mags.rds`, the
#'   depth of the MAGs with their GTDB taxonomy, and `singlem.rds`, the SingleM
#'   coverage of every lineage, both with the table of samples.
#' * `README.txt` lists every file with its rows, columns and content, and
#'   explains the units.
#'
#' Tables the run does not have (CAT without `--cat_db`, MAGs with
#' `--workflow contigs`...) are left out. The contig depth tables are not
#' written: they are large and are already plain tables in the pipeline
#' output.
#'
#' With several runs ([read_pipeline()] with several directories) the tables
#' keep their `study` column, and the samples and bins of the matrices are
#' named `<study>__<name>`.
#'
#' @param run A list returned by [read_pipeline()], or the results directory
#'   of a run (read with `depth = FALSE`).
#' @param output_dir Directory to write into; created if needed. Files with
#'   the same names are overwritten.
#' @param excel If `TRUE`, also write `magreport_tables.xlsx` (needs writexl).
#' @param phyloseq If `TRUE`, also write phyloseq objects (needs phyloseq).
#'
#' @return Invisibly, a data frame with one row per file written: `file`,
#'   `rows`, `columns` and `description`, the same as in `README.txt`.
#' @export
#' @seealso [mag_report()], which calls it with `tables = TRUE`.
#'
#' @examples
#' run <- read_pipeline(system.file("extdata", "example_contigs", package = "magreport"))
#' out <- file.path(tempdir(), "tables")
#' files <- export_tables(run, out, excel = FALSE, phyloseq = FALSE)
#' head(files)
#' head(utils::read.delim(file.path(out, "matrices", "singlem_phylum_relative.tsv"),
#'                        check.names = FALSE))
export_tables <- function(run, output_dir = "magreport_tables", excel = TRUE, phyloseq = TRUE) {
  if (is.character(run)) {
    run <- read_pipeline(run, depth = FALSE)
  }
  if (!is.list(run) || !all(c("mags", "mapping", "singlem") %in% names(run))) {
    stop("run must be a list returned by read_pipeline() or a results directory", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  contents <- new_contents()
  notes <- character()
  add <- function(x, file, description) {
    if (is.null(x) || nrow(x) == 0) {
      return(invisible(FALSE))
    }
    path <- file.path(output_dir, file)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write_table_tsv(x, path)
    contents[nrow(contents) + 1, ] <<- list(file, nrow(x), ncol(x), description)
    tables[[file]] <<- x
    invisible(TRUE)
  }
  tables <- list()

  # Summaries -----------------------------------------------------------------
  samples <- summary_samples(run)
  add(samples, "samples.tsv",
      "One row per sample: read pairs, mapping rate, gene catalog TPM and annotated fractions, SingleM coverage and MAGs with reads")
  add(summary_assemblies(run), "assemblies.tsv",
      "One row per assembly: contigs, length and N50 of each contig set, bins of each binner and quality of the final bins")
  bins <- summary_bins(run)
  add(bins, "bins.tsv",
      "One row per final bin (DAS Tool): size, CheckM2 quality, GTDB-Tk taxonomy, GUNC and annotation; is_representative marks the dereplicated MAGs")
  mags <- summary_mags(run, bins)
  add(mags, "mags.tsv",
      "One row per MAG (dereplicated representatives, or bins of medium or high quality): quality, taxonomy, abundance across samples and functions")

  # Methods -------------------------------------------------------------------
  methods <- run_methods(run)
  for (name in names(methods)) {
    add(methods[[name]], file.path("methods", paste0(name, ".tsv")),
        switch(name,
               run = "Pipeline version, run name, dates and command",
               databases = "Databases used, with their path and modification date",
               software = "Version of each tool, as reported during the run",
               params = "Parameters of the run"))
  }
  if (length(methods) == 0) {
    notes <- c(notes, "No manifest.json or run_info.txt found: methods/ was not written.")
  }

  # Matrices ------------------------------------------------------------------
  ab <- run$abundance
  if (!is.null(ab)) {
    add(wide_table(key_of(ab, "bin"), key_of(ab, "sample"), ab$mean_depth, "bin"),
        "matrices/mag_depth.tsv",
        "MAG x sample: mean read depth (coverage) of the MAG in the sample")
    add(wide_table(key_of(ab, "bin"), key_of(ab, "sample"), as.numeric(ab$rel_abundance), "bin"),
        "matrices/mag_relative.tsv",
        "MAG x sample: % of the summed depth of all the MAGs in the sample (columns add up to 100)")
    if ("reads_pct" %in% names(ab)) {
      add(wide_table(key_of(ab, "bin"), key_of(ab, "sample"), ab$reads_pct, "bin"),
          "matrices/mag_reads_pct.tsv",
          "MAG x sample: % of the reads of the sample that map to the MAG (CoverM relative abundance; the rest is unmapped)")
    }
    if ("covered_fraction" %in% names(ab)) {
      add(wide_table(key_of(ab, "bin"), key_of(ab, "sample"), ab$covered_fraction, "bin"),
          "matrices/mag_covered_fraction.tsv",
          "MAG x sample: fraction of the MAG covered by at least one read (0-1); low values mean a relative, not the MAG itself")
    }
    present <- tapply(ab$mean_depth > 0, key_of(ab, "bin"), sum)
    if (length(present) > 0 && all(present <= 1)) {
      notes <- c(notes, paste(
        "Every MAG has reads in one sample only: each sample was mapped to its own assembly",
        "(per_sample without --map_all_samples). The MAG matrices say in which sample each MAG",
        "was recovered and how abundant it is there, but they cannot compare a MAG between",
        "samples; for that use the CAT and SingleM profiles, or map every sample to the MAGs."))
    }
  }

  keep_mags <- if (is.null(mags)) NULL else key_of(mags, "bin")
  kos <- run$kos
  if (!is.null(kos)) {
    if (!is.null(keep_mags)) kos <- kos[key_of(kos, "bin") %in% keep_mags, , drop = FALSE]
    kos <- kos[!duplicated(data.frame(key_of(kos, "bin"), kos$ko)), , drop = FALSE]
    add(wide_table(kos$ko, key_of(kos, "bin"), kos$n_genes, "ko"),
        "matrices/mag_ko_genes.tsv",
        "KO x MAG: number of genes of the MAG annotated with the KEGG orthologue (eggNOG-mapper)")
  }
  hc <- run$hydrocarbon
  if (!is.null(hc)) {
    if (!is.null(keep_mags)) hc <- hc[key_of(hc, "bin") %in% keep_mags, , drop = FALSE]
    add(wide_table(paste(hc$database, hc$marker, sep = ":"), key_of(hc, "bin"), hc$n_genes, "marker"),
        "matrices/mag_hydrocarbon_genes.tsv",
        "Marker x MAG: genes of the MAG with each CANT-HYD or HADEG marker (database:marker)")
  }

  profiles <- list(cat = run$contig_taxonomy, singlem = run$singlem)
  for (source in names(profiles)) {
    x <- profiles[[source]]
    if (is.null(x)) next
    label <- if (source == "cat") "CAT (gene catalog)" else "SingleM (reads)"
    raw <- if (source == "cat") "tpm" else "coverage"
    for (rank in intersect(taxonomic_ranks(), unique(x$rank))) {
      r <- x[x$rank == rank, , drop = FALSE]
      if (!all(is.na(r[[raw]]))) {
        add(wide_table(r$taxon, key_of(r, "sample"), r[[raw]], "taxon"),
            sprintf("matrices/%s_%s_%s.tsv", source, rank, raw),
            sprintf("%s x sample, %s: %s", rank, label,
                    if (raw == "tpm") "summed TPM of the catalog genes of the taxon"
                    else "SingleM coverage of the taxon (reads per marker gene)"))
      }
      add(wide_table(r$taxon, key_of(r, "sample"), r$rel_abundance, "taxon"),
          sprintf("matrices/%s_%s_relative.tsv", source, rank),
          sprintf("%s x sample, %s: %% of the sample at that rank, Unclassified included (columns add up to 100)",
                  rank, label))
    }
  }

  features <- list(ko = run$contig_kos, canthyd = run$contig_canthyd, hadeg = run$contig_hadeg)
  feature_label <- c(ko = "KEGG orthologue (eggNOG-mapper)", canthyd = "CANT-HYD marker",
                     hadeg = "HADEG marker")
  for (name in names(features)) {
    x <- features[[name]]
    if (is.null(x)) next
    id <- if (name == "ko") "ko" else "marker"
    add(wide_table(x$feature, key_of(x, "sample"), x$tpm, id),
        sprintf("matrices/%s_tpm.tsv", name),
        sprintf("%s x sample, gene catalog: summed TPM of the genes with the %s", id, feature_label[[name]]))
    add(wide_table(x$feature, key_of(x, "sample"), x$rel_abundance, id),
        sprintf("matrices/%s_relative.tsv", name),
        sprintf("%s x sample, gene catalog: %% of the TPM of the sample in genes with the %s", id,
                feature_label[[name]]))
  }

  # Long tables ---------------------------------------------------------------
  long <- c(
    mapping = "Read mapping of each sample (Bowtie2)",
    assembly = "Contig statistics of each assembly and contig set",
    bins = "Every bin of every binner (CONCOCT, MetaBAT2, DAS Tool): size and GC",
    mags = "Final bins: CheckM2 quality and GTDB-Tk taxonomy, as read by read_mag_quality()",
    mimag = "Dereplicated MAGs: MIMAG quality with rRNA and tRNA genes, GUNC and GTDB-Tk",
    abundance = "MAG, sample, mean_depth and rel_abundance (% of the summed depth of the MAGs)",
    kos = "MAG, KEGG orthologue and number of genes",
    annotation = "Genes, annotated genes and KOs of each MAG",
    hydrocarbon = "MAG, database (CANT-HYD, HADEG), marker, category, pathway and number of genes",
    gunc = "GUNC chimera report of every bin set",
    catalog = "Size of the gene catalog",
    profile_stats = "Total TPM of each sample and fraction annotated by each source",
    contig_taxonomy = "CAT profile of the gene catalog: rank, taxon, sample, n_genes, tpm and rel_abundance",
    contig_kos = "KEGG profile of the gene catalog: ko, sample, n_genes, tpm and rel_abundance",
    contig_canthyd = "CANT-HYD profile of the gene catalog: marker, sample, n_genes, tpm and rel_abundance",
    contig_hadeg = "HADEG profile of the gene catalog: marker, sample, n_genes, tpm and rel_abundance",
    clusters = "Dereplication (galah): representative MAG and each member of its cluster",
    singlem = "SingleM profile of the reads: rank, taxon, sample, coverage and rel_abundance"
  )
  for (name in names(long)) {
    x <- run[[name]]
    if (!is.null(x) && "genome" %in% names(x) && "bin" %in% names(x) && identical(x$genome, x$bin)) {
      x$genome <- NULL
    }
    add(x, file.path("long", paste0(name, ".tsv")), long[[name]])
  }
  recovery <- mag_recovery_table(run)
  add(recovery, "long/mag_recovery.tsv",
      "Per sample and rank, % of the SingleM profile in taxa with at least one MAG of medium or high quality (mag_recovery())")

  # phyloseq ------------------------------------------------------------------
  if (isTRUE(phyloseq)) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      notes <- c(notes, "phyloseq is not installed: no phyloseq objects were written.")
      warning("Package 'phyloseq' is not installed: no phyloseq objects written", call. = FALSE)
    } else {
      ps <- list(
        mags = phyloseq_mags(run, mags, samples),
        singlem = phyloseq_singlem(run, samples)
      )
      for (name in names(ps)) {
        if (is.null(ps[[name]])) next
        file <- file.path("phyloseq", paste0(name, ".rds"))
        dir.create(file.path(output_dir, "phyloseq"), showWarnings = FALSE)
        saveRDS(ps[[name]], file.path(output_dir, file))
        contents[nrow(contents) + 1, ] <- list(
          file, phyloseq::ntaxa(ps[[name]]), phyloseq::nsamples(ps[[name]]),
          if (name == "mags") "phyloseq object: MAG depth (otu_table), GTDB-Tk taxonomy (tax_table) and samples.tsv (sample_data); rows = MAGs, columns = samples"
          else "phyloseq object: SingleM coverage of every lineage (otu_table), its taxonomy (tax_table) and samples.tsv (sample_data); rows = lineages, columns = samples"
        )
      }
      if (is.null(ps$singlem) && !is.null(run$singlem)) {
        notes <- c(notes, "SingleM was read from the tables per rank, without full lineages: phyloseq/singlem.rds was not written.")
      }
    }
  }
  if (!is.null(run$contig_taxonomy)) {
    notes <- c(notes, paste(
      "The CAT profile has one table per rank without the lineage of each taxon, so it has no phyloseq object:",
      "use matrices/cat_<rank>_*.tsv."))
  }

  # Excel ---------------------------------------------------------------------
  if (isTRUE(excel)) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      notes <- c(notes, "writexl is not installed: magreport_tables.xlsx was not written.")
      warning("Package 'writexl' is not installed: no Excel file written", call. = FALSE)
    } else {
      sheets_from <- contents$file[!startsWith(contents$file, "long/") &
                                     !startsWith(contents$file, "phyloseq/")]
      sheets <- list()
      for (file in sheets_from) {
        x <- tables[[file]]
        if (nrow(x) >= 1048576 || ncol(x) > 16384) {
          notes <- c(notes, sprintf("%s is larger than an Excel sheet: it is only in the TSV.", file))
          next
        }
        sheet <- sub("\\.tsv$", "", basename(file))
        if (startsWith(file, "methods/")) sheet <- paste0("methods_", sheet)
        sheets[[substr(sheet, 1, 31)]] <- x
      }
      xlsx <- "magreport_tables.xlsx"
      readme_sheet <- rbind(contents, list(xlsx, NA, NA, sprintf(
        "The summaries, methods and matrices, one sheet each (%d sheets, README first)",
        length(sheets) + 1)))
      writexl::write_xlsx(c(list(README = readme_sheet), sheets), file.path(output_dir, xlsx))
      contents <- readme_sheet
    }
  }

  write_readme(contents, notes, run, file.path(output_dir, "README.txt"))
  contents <- rbind(contents, list("README.txt", NA, NA, "This list of files"))
  invisible(contents)
}

# Helpers -------------------------------------------------------------------

new_contents <- function() {
  data.frame(file = character(), rows = integer(), columns = integer(),
             description = character(), stringsAsFactors = FALSE)
}

write_table_tsv <- function(x, path) {
  x[] <- lapply(x, function(column) {
    if (is.factor(column)) as.character(column) else as.vector(column)
  })
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     na = "NA", fileEncoding = "UTF-8")
}

# The name of a sample or bin, prefixed with its study when there are several
key_of <- function(x, column) {
  value <- as.character(x[[column]])
  if ("study" %in% names(x)) paste(x$study, value, sep = "__") else value
}

# A feature x sample table: first column `id`, then one column per sample
# (sorted), missing combinations as 0, rows by decreasing total
wide_table <- function(rows, columns, values, id) {
  if (length(rows) == 0) {
    return(NULL)
  }
  values <- as.numeric(values)
  m <- tapply(values, list(rows, columns), sum)
  m[is.na(m)] <- 0
  m <- m[order(-rowSums(m), rownames(m)), sort(colnames(m)), drop = FALSE]
  out <- data.frame(rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  names(out)[1] <- id
  rownames(out) <- NULL
  out
}

# Look up values of `from` for the rows of `to`, matching by study and `by`
lookup <- function(to, from, by, column) {
  if (is.null(from) || !column %in% names(from)) {
    return(rep(NA, nrow(to)))
  }
  from[[column]][match(key_of(to, by), key_of(from, by))]
}

# The ids of a table: its column `column`, and `study` when there are several
id_frame <- function(x, column) {
  out <- data.frame(as.character(x[[column]]), stringsAsFactors = FALSE)
  names(out) <- column
  if ("study" %in% names(x)) cbind(data.frame(study = x$study, stringsAsFactors = FALSE), out) else out
}

summary_samples <- function(run) {
  sources <- Filter(Negate(is.null), list(run$mapping, run$profile_stats, run$singlem,
                                          run$contig_taxonomy, run$abundance))
  if (length(sources) == 0) {
    return(NULL)
  }
  ids <- do.call(rbind, lapply(sources, function(x) {
    id_frame(x, "sample")
  }))
  out <- ids[!duplicated(ids), , drop = FALSE]
  out <- out[order(key_of(out, "sample")), , drop = FALSE]

  for (x in list(run$mapping, run$profile_stats)) {
    if (is.null(x)) next
    for (column in setdiff(names(x), c("study", "sample"))) {
      out[[column]] <- lookup(out, x, "sample", column)
    }
  }

  sm <- run$singlem
  if (!is.null(sm) && !all(is.na(sm$coverage))) {
    d <- sm[sm$rank == "domain", , drop = FALSE]
    total <- tapply(d$coverage, key_of(d, "sample"), sum)
    out$singlem_coverage <- as.numeric(total[key_of(out, "sample")])
  }

  ab <- run$abundance
  if (!is.null(ab)) {
    present <- tapply(ab$mean_depth > 0, key_of(ab, "sample"), sum)
    out$n_mags_with_reads <- as.integer(present[key_of(out, "sample")])
    if ("reads_pct" %in% names(ab)) {
      in_mags <- tapply(as.numeric(ab$reads_pct), key_of(ab, "sample"), sum)
      out$reads_in_mags_pct <- as.numeric(in_mags[key_of(out, "sample")])
    }
  }
  rownames(out) <- NULL
  out
}

summary_assemblies <- function(run) {
  sources <- Filter(Negate(is.null), list(run$assembly, run$bins, run$mags, run$mimag))
  sources <- Filter(function(x) "assembly" %in% names(x), sources)
  if (length(sources) == 0) {
    return(NULL)
  }
  ids <- do.call(rbind, lapply(sources, function(x) {
    id_frame(x, "assembly")
  }))
  out <- ids[!duplicated(ids) & !is.na(ids$assembly), , drop = FALSE]
  out <- out[order(key_of(out, "assembly")), , drop = FALSE]
  key <- key_of(out, "assembly")

  a <- run$assembly
  if (!is.null(a)) {
    a$set <- gsub("[^A-Za-z0-9]+", "", sub(">=", "min", a$contig_set))
    for (set in unique(a$set)) {
      s <- a[a$set == set, , drop = FALSE]
      for (column in intersect(c("n_contigs", "total_length", "n50", "max_length"), names(s))) {
        out[[paste(column, set, sep = "_")]] <- s[[column]][match(key, key_of(s, "assembly"))]
      }
    }
  }

  count_by <- function(x, group) {
    counts <- table(factor(key_of(x, "assembly"), levels = unique(key)), group)
    counts[key, , drop = FALSE]
  }
  if (!is.null(run$bins) && "binner" %in% names(run$bins)) {
    counts <- count_by(run$bins, run$bins$binner)
    for (binner in colnames(counts)) {
      out[[paste0("n_bins_", binner)]] <- as.integer(counts[, binner])
    }
  }
  if (!is.null(run$mags) && "quality" %in% names(run$mags)) {
    counts <- count_by(run$mags, factor(run$mags$quality, levels = quality_levels()))
    for (level in quality_levels()) {
      out[[paste0("n_", level)]] <- as.integer(counts[, level])
    }
  }
  if (!is.null(run$mimag) && "assembly" %in% names(run$mimag)) {
    out$n_representatives <- as.integer(count_by(run$mimag, rep("n", nrow(run$mimag)))[, 1])
  }
  rownames(out) <- NULL
  out
}

summary_bins <- function(run) {
  out <- if (!is.null(run$mags)) run$mags else run$bins
  if (is.null(out)) {
    return(NULL)
  }
  g <- run$gunc
  if (!is.null(g) && "bin" %in% names(g)) {
    g <- g[!duplicated(key_of(g, "bin")), , drop = FALSE]
    out$gunc_pass <- lookup(out, g, "bin", "pass")
    out$gunc_clade_separation_score <- lookup(out, g, "bin", "clade_separation_score")
    out$gunc_contamination_portion <- lookup(out, g, "bin", "contamination_portion")
    out$gunc_taxonomic_level <- lookup(out, g, "bin", "taxonomic_level")
  }
  ann <- run$annotation
  if (!is.null(ann)) {
    for (column in intersect(c("n_genes", "n_annotated", "n_kos", "annotated_pct"), names(ann))) {
      out[[column]] <- lookup(out, ann, "bin", column)
    }
  }
  if (!is.null(run$mimag)) {
    out$is_representative <- key_of(out, "bin") %in% key_of(run$mimag, "bin")
  } else if (!is.null(run$clusters)) {
    out$is_representative <- key_of(out, "bin") %in% representative_keys(run$clusters)
  }
  rownames(out) <- NULL
  out
}

# The dereplicated representatives, named as key_of() names bins
representative_keys <- function(clusters) {
  x <- clusters
  x$bin <- x$representative
  unique(key_of(x, "bin"))
}

summary_mags <- function(run, bins) {
  out <- run$mimag
  if (is.null(out)) {
    if (is.null(bins) || !"quality" %in% names(bins)) {
      return(NULL)
    }
    out <- if ("is_representative" %in% names(bins)) {
      bins[bins$is_representative, , drop = FALSE]
    } else {
      bins[bins$quality %in% c("high", "medium"), , drop = FALSE]
    }
    out$is_representative <- NULL
    if (nrow(out) == 0) {
      return(NULL)
    }
  } else {
    if ("genome" %in% names(out) && identical(out$genome, out$bin)) {
      out$genome <- NULL
    }
    first <- intersect(c("study", "bin"), names(out))
    out <- out[c(first, setdiff(names(out), first))]
    g <- run$gunc
    if (!is.null(g) && "bin" %in% names(g)) {
      g <- g[!duplicated(key_of(g, "bin")), , drop = FALSE]
      out$gunc_clade_separation_score <- lookup(out, g, "bin", "clade_separation_score")
    }
    ann <- run$annotation
    if (!is.null(ann)) {
      for (column in intersect(c("n_genes", "n_kos", "annotated_pct"), names(ann))) {
        out[[column]] <- lookup(out, ann, "bin", column)
      }
    }
  }
  key <- key_of(out, "bin")

  cl <- run$clusters
  if (!is.null(cl)) {
    cl$bin <- cl$representative
    members <- table(factor(key_of(cl, "bin"), levels = key))
    out$n_cluster_members <- as.integer(members)
  }

  ab <- run$abundance
  if (!is.null(ab)) {
    ab <- ab[key_of(ab, "bin") %in% key, , drop = FALSE]
    bk <- factor(key_of(ab, "bin"), levels = key)
    rel <- as.numeric(ab$rel_abundance)
    out$abundance_mean_pct <- as.numeric(tapply(rel, bk, mean))
    out$abundance_max_pct <- as.numeric(tapply(rel, bk, max))
    out$depth_mean <- as.numeric(tapply(ab$mean_depth, bk, mean))
    out$n_samples_with_reads <- as.integer(tapply(ab$mean_depth > 0, bk, sum))
    if ("reads_pct" %in% names(ab)) {
      out$reads_mean_pct <- as.numeric(tapply(as.numeric(ab$reads_pct), bk, mean))
      out$reads_max_pct <- as.numeric(tapply(as.numeric(ab$reads_pct), bk, max))
    }
  }

  if (!is.null(run$kos) && !"n_kos" %in% names(out)) {
    k <- run$kos[key_of(run$kos, "bin") %in% key, , drop = FALSE]
    n <- tapply(k$ko, factor(key_of(k, "bin"), levels = key), function(v) length(unique(v)))
    out$n_kos <- as.integer(ifelse(is.na(n), 0L, n))
  }

  hc <- run$hydrocarbon
  if (!is.null(hc)) {
    hc <- hc[key_of(hc, "bin") %in% key, , drop = FALSE]
    for (db in unique(hc$database)) {
      h <- hc[hc$database == db, , drop = FALSE]
      bk <- factor(key_of(h, "bin"), levels = key)
      markers <- tapply(h$marker, bk, function(v) paste(sort(unique(v)), collapse = ","))
      name <- tolower(gsub("[^A-Za-z0-9]+", "", db))
      out[[paste0("n_", name, "_markers")]] <- as.integer(tapply(h$marker, bk, function(v) length(unique(v))))
      out[[paste0("n_", name, "_markers")]][is.na(out[[paste0("n_", name, "_markers")]])] <- 0L
      out[[paste0(name, "_markers")]] <- as.character(markers)
    }
  }
  rownames(out) <- NULL
  out
}

run_methods <- function(run) {
  info <- attr(run, "run_info")
  if (is.null(info)) {
    return(list())
  }
  studies <- attr(run, "studies")
  if (is.null(studies)) {
    return(tryCatch(describe_run(info), error = function(e) list()))
  }
  per_study <- lapply(studies, function(study) {
    d <- tryCatch(describe_run(info[[study]]), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    lapply(d, function(x) {
      if (nrow(x) == 0) x else cbind(data.frame(study = study, stringsAsFactors = FALSE), x)
    })
  })
  per_study <- Filter(Negate(is.null), per_study)
  if (length(per_study) == 0) {
    return(list())
  }
  out <- lapply(names(per_study[[1]]), function(name) {
    do.call(rbind, lapply(per_study, function(d) {
      x <- d[[name]]
      if (!"study" %in% names(x)) NULL else x
    }))
  })
  names(out) <- names(per_study[[1]])
  Filter(Negate(is.null), out)
}

mag_recovery_table <- function(run) {
  mags <- if (!is.null(run$mimag)) run$mimag else run$mags
  if (is.null(run$singlem) || is.null(mags) || !"phylum" %in% names(mags)) {
    return(NULL)
  }
  tryCatch(mag_recovery(run$singlem, mags), error = function(e) NULL)
}

# sample_data of phyloseq: samples.tsv with the sample names as row names
phyloseq_samples <- function(samples, names) {
  if (is.null(samples)) {
    samples <- data.frame(sample = names, stringsAsFactors = FALSE)
  }
  rownames(samples) <- key_of(samples, "sample")
  phyloseq::sample_data(samples[intersect(names, rownames(samples)), , drop = FALSE])
}

phyloseq_mags <- function(run, mags, samples) {
  ab <- run$abundance
  if (is.null(ab) || is.null(mags)) {
    return(NULL)
  }
  depth <- wide_table(key_of(ab, "bin"), key_of(ab, "sample"), ab$mean_depth, "bin")
  otu <- as.matrix(depth[-1])
  rownames(otu) <- depth$bin

  ranks <- taxonomic_ranks()
  tax <- matrix(NA_character_, nrow(otu), length(ranks), dimnames = list(rownames(otu), ranks))
  hit <- match(rownames(otu), key_of(mags, "bin"))
  for (rank in intersect(ranks, names(mags))) {
    tax[, rank] <- as.character(mags[[rank]][hit])
  }
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax),
                     phyloseq_samples(samples, colnames(otu)))
}

# SingleM coverage of every full lineage, from the profiles of each sample
phyloseq_singlem <- function(run, samples) {
  files <- attr(run, "files")
  studies <- attr(run, "studies")
  paths <- if (is.null(studies)) list(files$singlem) else lapply(files, `[[`, "singlem")
  prefix <- if (is.null(studies)) "" else studies

  rows <- lapply(seq_along(paths), function(i) {
    p <- paths[[i]]
    study <- prefix[i]
    if (is.null(p)) return(NULL)
    tables <- lapply(p, read_pipeline_table, check.names = FALSE, quote = "")
    if (!all(vapply(tables, function(t) all(c("sample", "coverage", "taxonomy") %in% names(t)),
                    logical(1)))) {
      return(NULL)
    }
    x <- do.call(rbind, lapply(tables, `[`, c("sample", "coverage", "taxonomy")))
    x$coverage <- as.numeric(x$coverage)
    x <- x[!is.na(x$coverage), , drop = FALSE]
    if (study != "") x$sample <- paste(study, x$sample, sep = "__")
    x
  })
  x <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(x) || nrow(x) == 0) {
    return(NULL)
  }
  x$taxonomy <- sub("^Root;\\s*", "", trimws(x$taxonomy))
  wide <- wide_table(x$taxonomy, x$sample, x$coverage, "taxonomy")
  otu <- as.matrix(wide[-1])
  ids <- sprintf("singlem_%0*d", nchar(nrow(otu)), seq_len(nrow(otu)))
  rownames(otu) <- ids
  tax <- split_lineage(wide$taxonomy)
  rownames(tax) <- ids
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax),
                     phyloseq_samples(samples, colnames(otu)))
}

write_readme <- function(contents, notes, run, path) {
  studies <- attr(run, "studies")
  info <- attr(run, "run_info")
  version <- function(i) if (is.null(i$pipeline$version)) "unknown" else i$pipeline$version
  pipeline <- if (is.null(studies)) {
    if (is.null(info)) "unknown" else version(info)
  } else {
    paste(sprintf("%s: %s", studies, vapply(studies, function(s) version(info[[s]]), "")),
          collapse = "; ")
  }

  size <- ifelse(is.na(contents$rows), "",
                 sprintf("%9s x %-5s", format(contents$rows, big.mark = ","), contents$columns))
  files <- sprintf("%-42s %-17s %s", contents$file, size, contents$description)

  lines <- c(
    "Tables of a mag-pipeline run, written by magreport::export_tables()",
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M"),
           "   magreport ", as.character(utils::packageVersion("magreport")),
           "   mag-pipeline ", pipeline),
    "",
    "All tables are tab-separated, UTF-8, with a header row; missing values are NA.",
    "",
    "HOW TO READ THEM",
    "  R:      read.delim(\"matrices/mag_relative.tsv\", check.names = FALSE)",
    "  Python: pandas.read_csv(\"matrices/mag_relative.tsv\", sep=\"\\t\", index_col=0)",
    "  Excel:  magreport_tables.xlsx has the summaries, methods and matrices",
    "  phyloseq: readRDS(\"phyloseq/mags.rds\")",
    "",
    "FOLDERS",
    "  *.tsv      summaries, one row per sample, assembly, bin or MAG",
    "  matrices/  one row per MAG, taxon, KO or marker and one column per sample (or MAG),",
    "             ready for vegan, phyloseq, DESeq2, scikit-bio, Excel...",
    "  long/      the same data with one row per feature and sample (tidy), for ggplot2 or a database",
    "  methods/   pipeline version, databases, tool versions and parameters, for the methods of a paper",
    "  phyloseq/  phyloseq objects (R)",
    "",
    "UNITS: every abundance is given twice",
    "  mag_depth            mean read depth of the MAG in the sample: the depth of its contigs",
    "                       weighted by their length. A sample only has depth in the MAGs of the",
    "                       assemblies its reads were mapped to (all of them with coassembly or",
    "                       --map_all_samples; only its own assembly otherwise).",
    "  mag_relative         % of the summed depth of all the MAGs in the sample. It is the share",
    "                       among the MAGs, not the % of the reads of the sample.",
    "  mag_reads_pct        (only when every sample was mapped to the MAGs with CoverM) % of the",
    "                       reads of the sample that map to the MAG; reads_in_mags_pct in",
    "                       samples.tsv is their sum, the share of the sample the MAGs explain.",
    paste0("  Source of the MAG abundance: ", abundance_source(run), "."),
    "  cat_*_tpm, ko_tpm... summed TPM of the genes of the gene catalog (all contigs, binned or not)",
    "                       with that taxon or function; TPM already corrects for gene length and",
    "                       sequencing depth (each sample adds up to 1,000,000 over all genes)",
    "  cat_*_relative       % of the TPM of the sample at that rank, Unclassified included (adds up to 100)",
    "  ko/canthyd/hadeg_relative  % of the TPM of all the genes of the sample. A gene with several",
    "                       KOs counts for each, so KOs do not add up to 100.",
    "  singlem_*_coverage   SingleM coverage (reads per single-copy marker gene) of the taxon",
    "  singlem_*_relative   % of the SingleM coverage of the sample at that rank, Unclassified included",
    "  Raw read counts per gene are in the gene abundance tables of the pipeline (08_contigs),",
    "  which are not copied here.",
    "",
    "Taxon names follow GTDB (prefixes d__, p__, c__, o__, f__, g__, s__). \"Unclassified\" holds",
    "what was not classified down to that rank.",
    if (!is.null(studies)) c("",
      "SEVERAL STUDIES: tables have a column study; in the matrices samples and bins are named",
      "<study>__<name>."),
    "",
    "FILES (rows x columns)",
    files
  )
  if (length(notes) > 0) {
    lines <- c(lines, "", "NOTES", paste("  -", notes))
  }
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}

# Which table the MAG abundance came from, for the README
abundance_source <- function(run) {
  files <- attr(run, "files")
  paths <- if (is.null(attr(run, "studies"))) list(files$abundance) else lapply(files, `[[`, "abundance")
  paths <- unlist(paths)
  names <- if (is.null(paths)) character() else unique(basename(paths))
  if (length(names) == 0) {
    return("none")
  }
  describe <- ifelse(grepl("^mag_abundance_all_samples", names),
                     "every sample mapped to the dereplicated MAGs",
                     "the mapping of the pipeline")
  paste(sprintf("%s (%s)", names, describe), collapse = "; ")
}
