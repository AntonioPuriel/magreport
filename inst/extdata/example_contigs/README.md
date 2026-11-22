# Example output of the contig branch

A small **simulated** run of the contig branch of mag-pipeline 2.0.0
(`--workflow contigs` or `both`, with `--cat_db`, `--canthyd_db`,
`--eggnog_db` and `--singlem_db`), made by `data-raw/example_contigs.R`. The
file names, layout and columns are those written by the pipeline; the numbers
are simulated and mean nothing biologically.

12 samples named `<time>-<site>-<condition>-<replicate>`: times T0, T1 and T2,
site BB, a low-impact (C) and a high-impact (O) condition, 2 replicates.
Hydrocarbon degraders (Alcanivorax, Marinobacter) bloom at T1 in O.

- `profiles/taxonomy_per_sample.tsv.gz`: TPM of the catalog genes of each taxon
  (CAT) in each sample, domain to species.
- `profiles/ko_per_sample.tsv.gz`, `profiles/canthyd_per_sample.tsv.gz`: TPM of
  KEGG orthologues and CANT-HYD markers.
- `profiles/profile_stats.tsv`: share of the TPM of each sample annotated by
  each source.
- `catalog/gene_catalog.stats.tsv`: size of the gene catalog.
- `singlem/<sample>.singlem.tsv`: SingleM profile of the reads of each sample.
