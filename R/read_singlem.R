#' Read the SingleM profile of the reads
#'
#' Reads the taxonomic profiles written by SingleM in `09_singlem/`
#' (mag-pipeline 2.0.0 or later, with `--singlem_db`). SingleM classifies reads
#' of single-copy marker genes, so its profile needs neither assembly nor
#' binning and sees the organisms that never made it into a contig.
#'
#' Three layouts are recognised:
#'
#' * a profile with the columns `sample`, `coverage` and `taxonomy`: the
#'   profile of one sample (`profiles/<sample>.singlem.tsv`) or of all of them
#'   (`singlem_profiles_with_extras.tsv`). Every rank is computed from it: the
#'   coverage of a taxon is that of its reads plus those of its descendants,
#'   and the coverage of reads classified only to a higher rank counts as
#'   `"Unclassified"`.
#' * a taxon-by-sample table of relative abundances at one rank
#'   (`singlem_<rank>.tsv`). Lineages that stop above that rank count as
#'   `"Unclassified"`.
#'
#' @param path Path to one or more SingleM tables of the same layout, for
#'   example every file of `09_singlem/profiles/`.
#' @param rank For a taxon-by-sample table, its rank. By default it is taken
#'   from the file name (`singlem_genus.tsv`) or else from the deepest GTDB
#'   prefix found.
#'
#' @return A data frame with the same columns as [read_taxonomy_profile()]:
#'   `rank`, `taxon`, `sample`, `coverage` (`NA` for a taxon-by-sample table)
#'   and `rel_abundance`, the percentage of the sample at that rank. Taxa keep
#'   their GTDB prefix (`p__Pseudomonadota`). The attribute `"source"` is
#'   `"SingleM"`.
#' @export
#' @seealso [plot_community()], [plot_profile_comparison()].
#'
#' @examples
#' dir <- system.file("extdata", "example_contigs", "singlem", package = "magreport")
#' singlem <- read_singlem(list.files(dir, full.names = TRUE))
#' head(singlem[singlem$rank == "phylum", ])
read_singlem <- function(path, rank = NULL) {
  if (length(path) == 0) {
    stop("No SingleM table given", call. = FALSE)
  }
  # Sample names such as T0-BB-C-1 are column names in the taxon-by-sample
  # tables: keep them as they are
  tables <- lapply(path, read_pipeline_table, check.names = FALSE, quote = "")
  first <- names(tables[[1]])

  if (all(c("sample", "coverage", "taxonomy") %in% first)) {
    x <- do.call(rbind, lapply(tables, function(t) {
      t[c("sample", "coverage", "taxonomy")]
    }))
    out <- singlem_from_coverage(x)
  } else if ("taxonomy" %in% first) {
    out <- do.call(rbind, lapply(seq_along(tables), function(i) {
      singlem_from_wide(tables[[i]], path[i], rank)
    }))
  } else {
    stop("Unexpected format in ", basename(path[1]),
         ": expected the column taxonomy", call. = FALSE)
  }

  out <- out[order(match(out$rank, taxonomic_ranks()), out$sample, -out$rel_abundance), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "source") <- "SingleM"
  out
}

# "Root; d__Bacteria; p__Bacillota" -> c(domain = "d__Bacteria", phylum = "p__Bacillota", ...)
split_lineage <- function(taxonomy) {
  prefixes <- rank_prefixes()
  parts <- strsplit(as.character(taxonomy), ";", fixed = TRUE)
  out <- t(vapply(parts, function(p) {
    p <- trimws(p)
    vapply(prefixes, function(pre) {
      hit <- p[startsWith(p, pre)]
      if (length(hit) == 0 || hit[1] == pre) NA_character_ else hit[1]
    }, character(1))
  }, character(length(prefixes))))
  colnames(out) <- names(prefixes)
  out
}

singlem_from_coverage <- function(x) {
  x$coverage <- as.numeric(x$coverage)
  x <- x[!is.na(x$coverage), , drop = FALSE]
  lineage <- split_lineage(x$taxonomy)

  out <- do.call(rbind, lapply(taxonomic_ranks(), function(r) {
    taxon <- lineage[, r]
    taxon[is.na(taxon)] <- "Unclassified"
    sums <- stats::aggregate(
      list(coverage = x$coverage),
      by = list(taxon = taxon, sample = x$sample),
      FUN = sum
    )
    sums$rank <- r
    sums
  }))

  out$rel_abundance <- share_within(out$coverage, paste(out$sample, out$rank, sep = "\r"))
  out[c("rank", "taxon", "sample", "coverage", "rel_abundance")]
}

singlem_from_wide <- function(x, path, rank) {
  samples <- setdiff(names(x), "taxonomy")
  if (length(samples) == 0) {
    stop("No sample columns in ", basename(path), call. = FALSE)
  }
  lineage <- split_lineage(x$taxonomy)

  if (is.null(rank)) {
    from_name <- sub("^singlem_([a-z]+)\\.tsv(\\.gz)?$", "\\1", basename(path))
    rank <- if (from_name %in% taxonomic_ranks()) {
      from_name
    } else {
      found <- taxonomic_ranks()[colSums(!is.na(lineage)) > 0]
      if (length(found) == 0) {
        stop("Could not tell the rank of ", basename(path), ": use rank", call. = FALSE)
      }
      found[length(found)]
    }
  }
  rank <- match.arg(rank, taxonomic_ranks())

  taxon <- lineage[, rank]
  taxon[is.na(taxon)] <- "Unclassified"

  out <- do.call(rbind, lapply(samples, function(s) {
    value <- as.numeric(x[[s]])
    value[is.na(value)] <- 0
    data.frame(taxon = taxon, sample = s, value = value, stringsAsFactors = FALSE)
  }))
  out <- stats::aggregate(value ~ taxon + sample, data = out, FUN = sum)

  # Fractions or percentages; what is missing from 100 % was not classified
  # down to this rank
  totals <- tapply(out$value, out$sample, sum)
  scale <- if (all(totals <= 1.5)) 100 else 1
  out$value <- out$value * scale
  totals <- totals * scale
  gap <- 100 - totals
  gap <- gap[gap > 0.5]
  if (length(gap) > 0) {
    out <- rbind(out, data.frame(taxon = "Unclassified", sample = names(gap),
                                 value = as.numeric(gap), stringsAsFactors = FALSE))
    out <- stats::aggregate(value ~ taxon + sample, data = out, FUN = sum)
  }

  data.frame(
    rank = rank,
    taxon = out$taxon,
    sample = out$sample,
    coverage = NA_real_,
    rel_abundance = share_within(out$value, out$sample),
    stringsAsFactors = FALSE
  )
}
