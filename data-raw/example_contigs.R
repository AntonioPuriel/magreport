# Small simulated run of the contig branch of mag-pipeline (v2.0.0), for the
# examples and tests of magreport. The layout and columns are those written by
# bin/contig_profiles.py, modules/gene_catalog.nf and modules/singlem.nf of
# mag-pipeline; the numbers are simulated.
#
# 12 samples: times T0, T1 and T2, a low-impact (C) and a high-impact (O)
# condition, 2 replicates. Hydrocarbon degraders (Alcanivorax, Marinobacter)
# bloom at T1 in O.
#
#   Rscript data-raw/example_contigs.R    (from the package root)

set.seed(20261122)
out <- file.path("inst", "extdata", "example_contigs")
dir.create(file.path(out, "profiles"), recursive = TRUE, showWarnings = FALSE)

samples <- expand.grid(replicate = 1:2, condition = c("C", "O"), time = c("T0", "T1", "T2"),
                       stringsAsFactors = FALSE)
samples$sample <- with(samples, paste(time, "BB", condition, replicate, sep = "-"))

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Lineages, some of them known only down to a higher rank
lineages <- c(
  "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Alcanivoracaceae;g__Alcanivorax;s__Alcanivorax borkumensis",
  "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Marinobacteraceae;g__Marinobacter",
  "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas",
  "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Sulfitobacter",
  "d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae",
  "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;g__Polaribacter",
  "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae",
  "d__Bacteria;p__Desulfobacterota;c__Desulfobulbia;o__Desulfobulbales;f__Desulfocapsaceae;g__Desulfocapsa",
  "d__Bacteria;p__Desulfobacterota;c__Desulfobacteria",
  "d__Bacteria;p__Actinomycetota;c__Acidimicrobiia;o__Acidimicrobiales;f__Ilumatobacteraceae;g__Ilumatobacter",
  "d__Bacteria;p__Chloroflexota;c__Anaerolineae;o__Anaerolineales;f__Anaerolineaceae",
  "d__Bacteria;p__Planctomycetota",
  "d__Bacteria",
  "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosopumilus",
  "d__Archaea;p__Asgardarchaeota",
  ""   # contigs that CAT could not classify at all
)
base <- c(3, 4, 3, 6, 5, 7, 6, 5, 4, 3, 5, 3, 8, 1.5, 0.5, 30)
bloom <- c(12, 5, 2, rep(1, 13))   # fold change at T1 in O

split <- function(l) {
  parts <- if (l == "") character() else strsplit(l, ";", fixed = TRUE)[[1]]
  c(parts, rep(NA_character_, length(ranks) - length(parts)))
}
lineage_table <- t(vapply(lineages, split, character(length(ranks))))
colnames(lineage_table) <- ranks
rownames(lineage_table) <- NULL

share <- sapply(seq_len(nrow(samples)), function(i) {
  fold <- if (samples$time[i] == "T1" && samples$condition[i] == "O") bloom else
    if (samples$time[i] == "T2" && samples$condition[i] == "O") sqrt(bloom) else 1
  v <- base * fold * exp(stats::rnorm(length(base), 0, 0.15))
  v / sum(v)
})
tpm <- 1e6 * share                       # lineages x samples
n_genes <- round(base * 9000 * stats::runif(length(base), 0.8, 1.2))

# taxonomy_per_sample.tsv: rank, taxon, sample, n_genes, tpm
taxonomy <- do.call(rbind, lapply(ranks, function(r) {
  taxon <- lineage_table[, r]
  taxon[is.na(taxon)] <- "Unclassified"
  do.call(rbind, lapply(sort(unique(taxon)), function(t) {
    rows <- taxon == t
    data.frame(rank = r, taxon = t, sample = samples$sample,
               n_genes = sum(n_genes[rows]), tpm = colSums(tpm[rows, , drop = FALSE]))
  }))
}))
taxonomy$tpm <- signif(taxonomy$tpm, 6)
gz <- gzfile(file.path(out, "profiles", "taxonomy_per_sample.tsv.gz"), "w")
utils::write.table(taxonomy, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

# canthyd_per_sample.tsv: marker, description, sample, n_genes, tpm. Carried
# mostly by Alcanivorax and Marinobacter, so they follow their bloom.
markers <- data.frame(
  marker = c("AlkB", "CYP153", "AlmA_GroupI", "NdoB", "NdoC", "AhyA"),
  description = c("alkane hydrolase", "alkane oxidizing cytochrome P450",
                  "flavin-binding alkane monooxygenase",
                  "benzene/toluene/naphtalene dioxygenase subunit alpha",
                  "benzene/toluene/naphtalene dioxygenase subunit beta",
                  "molybdopterin-family alkane C2 methylene hydroxylase"),
  per_degrader = c(60, 90, 15, 20, 20, 0),
  background = c(15, 30, 5, 25, 25, 8),
  stringsAsFactors = FALSE
)
degraders <- colSums(tpm[1:2, ])
canthyd <- do.call(rbind, lapply(seq_len(nrow(markers)), function(m) {
  v <- (markers$per_degrader[m] * degraders / 1e4 + markers$background[m]) *
    exp(stats::rnorm(ncol(tpm), 0, 0.1))
  data.frame(marker = markers$marker[m], description = markers$description[m],
             sample = samples$sample, n_genes = round(stats::runif(1, 20, 300)),
             tpm = signif(v, 6))
}))
utils::write.table(canthyd, gzfile(file.path(out, "profiles", "canthyd_per_sample.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# ko_per_sample.tsv: ko, sample, n_genes, tpm. The KOs of hydrocarbon_kos() and
# some housekeeping genes
kos <- c("K00496", "K00114", "K00128", "K05549", "K03381", "K00446", "K14579",
         "K07540", "K03781", "K04564", "K11180", "K00394", "K02588", "K00370",
         "K02946", "K02886", "K02906", "K03043", "K03046", "K01869")
alkane <- kos %in% c("K00496", "K00114", "K00128")
ko <- do.call(rbind, lapply(seq_along(kos), function(k) {
  level <- if (k > 14) 900 else 120
  v <- (level + if (alkane[k]) 4 * degraders / 1e4 else 0) * exp(stats::rnorm(ncol(tpm), 0, 0.15))
  data.frame(ko = kos[k], sample = samples$sample,
             n_genes = round(stats::runif(1, 100, 3000)), tpm = signif(v, 6))
}))
utils::write.table(ko, gzfile(file.path(out, "profiles", "ko_per_sample.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# profile_stats.tsv
unclassified_domain <- tpm[nrow(tpm), ] / 1e6
stats <- data.frame(
  sample = samples$sample,
  total_tpm = 1e6,
  fraction_eggnog_ko = round(stats::runif(ncol(tpm), 0.52, 0.61), 4),
  fraction_canthyd = round(colSums(matrix(canthyd$tpm, ncol = ncol(tpm), byrow = TRUE)) / 1e6, 6),
  fraction_cat_domain = round(1 - unclassified_domain, 4)
)
utils::write.table(stats, file.path(out, "profiles", "profile_stats.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# gene_catalog.stats.tsv
dir.create(file.path(out, "catalog"), showWarnings = FALSE)
writeLines(c("genes_called\tgenes_kept\tcatalog_genes\tsingletons",
             "1254310\t1198772\t702115\t512004"),
           file.path(out, "catalog", "gene_catalog.stats.tsv"))

# SingleM: one profile per sample (sample, coverage, taxonomy), as written by
# `singlem pipe -p`. SingleM sees more archaea and leaves less unclassified.
dir.create(file.path(out, "singlem"), showWarnings = FALSE)
singlem_base <- base
singlem_base[14:15] <- c(6, 2)
singlem_base[16] <- 0
for (i in seq_len(nrow(samples))) {
  fold <- if (samples$time[i] == "T1" && samples$condition[i] == "O") bloom else
    if (samples$time[i] == "T2" && samples$condition[i] == "O") sqrt(bloom) else 1
  cov <- singlem_base * fold * exp(stats::rnorm(length(base), 0, 0.2))
  taxonomy_string <- vapply(lineages, function(l) {
    paste(c("Root", if (l != "") strsplit(l, ";", fixed = TRUE)[[1]]), collapse = "; ")
  }, character(1))
  keep <- cov > 0
  profile <- data.frame(sample = samples$sample[i], coverage = round(cov[keep], 2),
                        taxonomy = taxonomy_string[keep])
  profile <- rbind(profile, data.frame(sample = samples$sample[i], coverage = 0.8,
                                       taxonomy = "Root"))
  utils::write.table(profile, file.path(out, "singlem", paste0(samples$sample[i], ".singlem.tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}
