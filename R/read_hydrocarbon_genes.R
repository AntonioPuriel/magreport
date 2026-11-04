#' Read the hydrocarbon degradation genes of each MAG
#'
#' Reads `07_function/hydrocarbon_per_mag.tsv`, written by mag-pipeline
#' (version 1.9.0 or later) when it is run with `--canthyd_db` and/or
#' `--hadeg_db`. Each row is a marker found in a genome, with the database that
#' found it (`CANT-HYD` or `HADEG`), its category and, for HADEG, its pathway.
#'
#' @param path Path to `hydrocarbon_per_mag.tsv`.
#'
#' @return A data frame with one row per genome, database and marker. The
#'   column `bin` (a copy of `genome`) is added to match the other tables.
#' @export
#' @seealso [hydrocarbon_groups()] to group the CANT-HYD markers by process,
#'   [plot_hydrocarbon_genes()].
#'
#' @examples
#' path <- system.file("extdata", "marine_honour", "hydrocarbon_per_mag.tsv",
#'                     package = "magreport")
#' genes <- read_hydrocarbon_genes(path)
#' table(genes$database)
read_hydrocarbon_genes <- function(path) {
  x <- read_pipeline_table(path, required = c("genome", "database", "marker"))

  if ("n_genes" %in% names(x)) {
    x$n_genes <- as.numeric(x$n_genes)
  } else {
    x$n_genes <- 1
  }

  x$bin <- x$genome
  x
}

#' Group CANT-HYD markers by the process they take part in
#'
#' Assigns every CANT-HYD marker found in a MAG to a process: aerobic
#' degradation of alkanes, of short-chain alkanes or of aromatics and PAHs, or
#' anaerobic degradation. The result can be passed as `groups` to
#' [plot_abundance_timeline()] to follow each group of degraders over time.
#'
#' A MAG with markers of several processes appears once in each of them.
#' Markers that do not match any rule are put in `"Other"`.
#'
#' @param genes A data frame returned by [read_hydrocarbon_genes()].
#' @param any If `TRUE` (the default), adds the group `"Any CANT-HYD marker"`
#'   with every MAG that carries at least one marker.
#'
#' @return A data frame with the columns `bin` and `group`.
#' @export
#'
#' @details
#' The rules follow the gene names of the CANT-HYD HMM set:
#'
#' * *Alkanes (aerobic)*: `AlkB`, `CYP153`, `AlmA`, `LadA`, `LadB`.
#' * *Short-chain alkanes (aerobic)*: propane and butane monooxygenases
#'   (`PrmA`, `PrmC`, `pBmo`).
#' * *Aromatics and PAHs (aerobic)*: ring-hydroxylating dioxygenases (`NdoB`,
#'   `NdoC`, `MAH`, non-`NdoB` types) and toluene/phenol monooxygenases
#'   (`Tmo`, `Tom`, `Dmp`).
#' * *Anaerobic*: fumarate addition (`AssA`, `BssA`, `NmsA`), alkane C2
#'   hydroxylation (`AhyA`) and anaerobic aromatic activation (`AbcA`,
#'   `EbdA`, `CmdA`).
#'
#' @references Khot V, et al. (2022). CANT-HYD: a curated database of
#'   phylogeny-derived hidden Markov models for annotation of marker genes
#'   involved in hydrocarbon degradation. *Frontiers in Microbiology* 12:764058.
#'
#' @examples
#' path <- system.file("extdata", "marine_honour", "hydrocarbon_per_mag.tsv",
#'                     package = "magreport")
#' groups <- hydrocarbon_groups(read_hydrocarbon_genes(path))
#' table(groups$group)
hydrocarbon_groups <- function(genes, any = TRUE) {
  stopifnot(is.data.frame(genes), all(c("bin", "database", "marker") %in% names(genes)))

  canthyd <- genes[genes$database == "CANT-HYD", , drop = FALSE]
  if (nrow(canthyd) == 0) {
    stop("No CANT-HYD markers in genes: was the pipeline run with --canthyd_db?", call. = FALSE)
  }

  rules <- c(
    "Alkanes (aerobic)"             = "^(AlkB|CYP153|AlmA|LadA|LadB)",
    "Short-chain alkanes (aerobic)" = "^(PrmA|PrmC|pBmo)",
    "Aromatics and PAHs (aerobic)"  = "^(NdoB|NdoC|non_NdoB|MAH|Tmo|Tom|Dmp)",
    "Anaerobic"                     = "^(AhyA|AssA|BssA|NmsA|AbcA|EbdA|CmdA)"
  )

  group <- rep("Other", nrow(canthyd))
  for (name in names(rules)) {
    group[group == "Other" & grepl(rules[[name]], canthyd$marker)] <- name
  }

  out <- unique(data.frame(bin = canthyd$bin, group = group, stringsAsFactors = FALSE))

  if (isTRUE(any)) {
    out <- rbind(
      data.frame(bin = unique(canthyd$bin), group = "Any CANT-HYD marker"),
      out
    )
  }

  levels_used <- c("Any CANT-HYD marker", names(rules), "Other")
  out$group <- factor(out$group, levels = levels_used[levels_used %in% out$group])
  out <- out[order(out$group, out$bin), ]
  rownames(out) <- NULL
  out
}
