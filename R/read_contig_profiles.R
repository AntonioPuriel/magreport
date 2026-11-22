#' Read the taxonomic profile of the gene catalog
#'
#' Reads `08_contigs/profiles/taxonomy_per_sample.tsv`, written by mag-pipeline
#' (version 2.0.0 or later) with `--workflow contigs` or `both` and `--cat_db`.
#' Every catalog gene takes the CAT classification of its contig, and the table
#' holds the summed TPM of the genes of each taxon in each sample, from domain
#' to species. Genes on contigs that CAT could not classify at a rank are
#' counted as `"Unclassified"`.
#'
#' Unlike [read_mag_abundance()], this profile covers the whole community: the
#' organisms that were assembled but never binned count too.
#'
#' @param path Path to `taxonomy_per_sample.tsv`.
#'
#' @return A data frame with one row per rank, taxon and sample, with the
#'   columns `rank`, `taxon`, `sample`, `n_genes` (catalog genes of that
#'   taxon), `tpm` and `rel_abundance`: the percentage of the TPM of the sample
#'   at that rank, `"Unclassified"` included, so that every rank of every
#'   sample adds up to 100. The attribute `"source"` is `"CAT"`.
#' @export
#' @seealso [read_singlem()] for the profile from the reads, with the same
#'   columns; [plot_community()].
#'
#' @examples
#' path <- system.file("extdata", "example_contigs", "profiles",
#'                     "taxonomy_per_sample.tsv.gz", package = "magreport")
#' cat_profile <- read_taxonomy_profile(path)
#' head(cat_profile[cat_profile$rank == "phylum", ])
read_taxonomy_profile <- function(path) {
  x <- read_pipeline_table(path, required = c("rank", "taxon", "sample", "tpm"), quote = "")
  x$tpm <- as.numeric(x$tpm)
  if ("n_genes" %in% names(x)) {
    x$n_genes <- as.numeric(x$n_genes)
  }
  x$taxon[is.na(x$taxon) | x$taxon == ""] <- "Unclassified"

  unknown <- setdiff(unique(x$rank), taxonomic_ranks())
  if (length(unknown) > 0) {
    stop("Unexpected rank(s) in ", basename(path), ": ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }

  x$rel_abundance <- share_within(x$tpm, paste(x$sample, x$rank, sep = "\r"))
  x <- x[order(match(x$rank, taxonomic_ranks()), x$sample, -x$tpm), , drop = FALSE]
  rownames(x) <- NULL
  attr(x, "source") <- "CAT"
  x
}

#' Read a functional profile of the gene catalog
#'
#' Reads the profiles of KEGG orthologues and hydrocarbon markers in
#' `08_contigs/profiles/` (mag-pipeline 2.0.0 or later, `--workflow contigs` or
#' `both`): `ko_per_sample.tsv` (eggNOG-mapper), `canthyd_per_sample.tsv` and
#' `hadeg_per_sample.tsv`. Each row is the summed TPM of the catalog genes
#' carrying a feature in a sample. A gene with several KOs counts towards each.
#'
#' @param path Path to one of the `*_per_sample.tsv` tables.
#' @param total_tpm Total TPM of each sample, the denominator of
#'   `rel_abundance`: a named numeric vector (names are samples) or the data
#'   frame returned by [read_profile_stats()]. By default 10^6, the TPM of
#'   every catalog gene together.
#'
#' @return A data frame with one row per feature and sample: the columns of
#'   the table (`ko`, or `marker` with its description, compound or pathway),
#'   `feature` (a copy of the first of them), `sample`, `n_genes`, `tpm` and
#'   `rel_abundance`, the percentage of the TPM of the sample.
#' @export
#' @seealso [plot_feature_heatmap()], [plot_abundance_timeline()].
#'
#' @examples
#' path <- system.file("extdata", "example_contigs", "profiles",
#'                     "canthyd_per_sample.tsv.gz", package = "magreport")
#' markers <- read_feature_profile(path)
#' head(markers)
read_feature_profile <- function(path, total_tpm = NULL) {
  x <- read_pipeline_table(path, required = c("sample", "tpm"), quote = "")
  x$tpm <- as.numeric(x$tpm)
  if ("n_genes" %in% names(x)) {
    x$n_genes <- as.numeric(x$n_genes)
  }

  features <- setdiff(names(x), c("sample", "n_genes", "tpm"))
  if (length(features) == 0) {
    stop("No feature column (ko, marker...) in ", basename(path), call. = FALSE)
  }
  x$feature <- as.character(x[[features[1]]])

  if (is.data.frame(total_tpm)) {
    if (!all(c("sample", "total_tpm") %in% names(total_tpm))) {
      stop("total_tpm must come from read_profile_stats()", call. = FALSE)
    }
    total_tpm <- stats::setNames(total_tpm$total_tpm, total_tpm$sample)
  }
  total <- if (is.null(total_tpm)) {
    rep(1e6, nrow(x))
  } else {
    missing <- setdiff(unique(x$sample), names(total_tpm))
    if (length(missing) > 0) {
      stop("total_tpm has no value for sample(s): ",
           paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
    }
    unname(total_tpm[x$sample])
  }
  x$rel_abundance <- ifelse(total > 0, 100 * x$tpm / total, 0)

  x <- x[c(features, "feature", setdiff(names(x), c(features, "feature")))]
  rownames(x) <- NULL
  x
}

#' Read how much of each sample the annotations explain
#'
#' Reads `08_contigs/profiles/profile_stats.tsv`: for every sample, the total
#' TPM of the gene catalog and the fraction of it carried by genes with a KEGG
#' orthologue, a CANT-HYD or HADEG marker, or a contig classified by CAT.
#'
#' @param path Path to `profile_stats.tsv`.
#'
#' @return A data frame with one row per sample and the columns `sample`,
#'   `total_tpm` and one `fraction_<source>` column per annotation that was
#'   run, between 0 and 1.
#' @export
#' @seealso [plot_annotation_fraction()].
#'
#' @examples
#' path <- system.file("extdata", "example_contigs", "profiles", "profile_stats.tsv",
#'                     package = "magreport")
#' read_profile_stats(path)
read_profile_stats <- function(path) {
  x <- read_pipeline_table(path, required = c("sample", "total_tpm"))
  for (column in c("total_tpm", grep("^fraction_", names(x), value = TRUE))) {
    x[[column]] <- as.numeric(x[[column]])
  }
  x
}

#' Read the size of the gene catalog
#'
#' Reads `08_contigs/catalog/gene_catalog.stats.tsv`: genes called by Prodigal
#' in every assembly, genes long enough to enter the catalog, non-redundant
#' catalog genes after clustering (95 % identity over 90 % of the shorter gene)
#' and how many of them are singletons.
#'
#' @param path Path to `gene_catalog.stats.tsv`.
#'
#' @return A one-row data frame with numeric columns.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "example_contigs", "catalog", "gene_catalog.stats.tsv",
#'                     package = "magreport")
#' read_catalog_stats(path)
read_catalog_stats <- function(path) {
  x <- read_pipeline_table(path, required = c("genes_called", "catalog_genes"))
  x[] <- lapply(x, as.numeric)
  x
}

# Ranks used by CAT, SingleM and GTDB-Tk, and their GTDB prefixes
taxonomic_ranks <- function() {
  c("domain", "phylum", "class", "order", "family", "genus", "species")
}

rank_prefixes <- function() {
  stats::setNames(c("d__", "p__", "c__", "o__", "f__", "g__", "s__"), taxonomic_ranks())
}

# Percentage of each value within its group
share_within <- function(value, group) {
  totals <- tapply(value, group, sum)
  total <- as.numeric(totals[group])
  ifelse(total > 0, 100 * value / total, 0)
}
