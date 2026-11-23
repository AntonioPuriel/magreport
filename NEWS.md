# magreport (development version)

* `read_manifest()` reads `manifest.json`, which mag-pipeline 2.1.0 and later
  write in the results directory: pipeline version, run, databases, parameters
  and the list of tables of the run. `read_pipeline()` uses it when it is
  there, so it no longer guesses file names, and keeps it in the attribute
  `"manifest"`. Without a manifest (mag-pipeline 2.0.0 or earlier) the tables
  are still found by name. magreport reads `schema_version` 1.x.
* `read_run_info()` reads `run_info.txt` of older runs into the same form, and
  `describe_run()` turns either into tables for the methods of a paper:
  versions, dates, databases with their modification date and parameters.
* `read_gunc()` reads the GUNC chimera reports (`--gunc_db`) and `plot_gunc()`
  draws the clade separation score of every bin against its CheckM2
  contamination. `read_pipeline()` returns them in `gunc`.
* `mag_report()` has two new sections: "Chimeras (GUNC)", with the bins that
  CheckM2 rates high or medium but GUNC finds chimeric, and "Methods", with
  the pipeline, databases and parameters of the run.
* `mag_recovery()` computes, per sample and rank, the share of the SingleM (or
  CAT) profile in taxa with at least one MAG, and `plot_mag_recovery()` draws
  it. `mag_report()` shows it in "Community profile".
* New example data: the GUNC reports of the 1161 bins of Marine Honour
  (`extdata/marine_honour/gunc`).
* New dependency: jsonlite.

# magreport 0.8.0

## The community without MAGs

mag-pipeline 2.0 profiles the whole community besides recovering MAGs: a gene
catalog of every contig with its abundance (TPM), taxonomy (CAT) and function
(eggNOG, CANT-HYD, HADEG) in each sample, and a SingleM profile of the reads.
magreport now reads and draws these tables.

* `read_taxonomy_profile()` reads `taxonomy_per_sample.tsv`, the CAT profile of
  the gene catalog, and adds the relative abundance of each taxon at each rank.
* `read_singlem()` reads the SingleM profiles (of each sample, of all samples,
  or the taxon-by-sample tables per rank) into the same columns, every rank
  computed from the coverage of the reads.
* `read_feature_profile()` reads the KEGG orthologue and hydrocarbon marker
  profiles (`ko_per_sample.tsv`, `canthyd_per_sample.tsv`,
  `hadeg_per_sample.tsv`); `read_profile_stats()` and `read_catalog_stats()`
  read the annotated share of each sample and the size of the catalog.
* `plot_community()` draws the composition of every sample at any rank, from
  CAT or SingleM. `plot_profile_comparison()` compares the two views taxon by
  taxon. `plot_annotation_fraction()` shows how much of each sample the
  annotations explain, and `plot_feature_heatmap()` the abundance of genes or
  markers, for example those of `hydrocarbon_kos()`.
* `plot_abundance_timeline()` also follows genes, markers and taxa over time,
  not only MAGs.
* `read_pipeline()` finds the new tables and `mag_report()` adds a
  "Community profile" section; it no longer needs MAGs, so a run with
  `--workflow contigs` gets a report too.
* New example data `example_contigs`: a small simulated run of the contig
  branch, in the format of mag-pipeline 2.0.

# magreport 0.7.0

* `plot_mag_tree(interactive = TRUE)` returns an interactive version of the
  tree (a ggiraph htmlwidget). Pointing at a MAG shows its name, taxonomy,
  MIMAG quality, CANT-HYD markers, highest abundance and response, and
  highlights it in every ring. Pointing at a coloured cell of a hydrocarbon ring
  highlights every MAG of that process. A search box above the tree highlights
  and lists the MAGs whose name, taxonomy, lineage or CANT-HYD markers contain
  the text typed (for example "Alcanivorax" or "AhyA"). ggiraph is suggested.

# magreport 0.6.0

* `plot_mag_tree()` draws the MAGs as a circular tree with iTOL-style rings:
  lineage, novelty against GTDB, MIMAG quality (with its thresholds), one ring
  per hydrocarbon degradation process, the response to a disturbance (log2 of
  the mean abundance between two conditions) and bars with the highest
  abundance of each MAG, on a square-root scale by default. Every ring after the
  lineage is optional. It needs ggtree (Bioconductor), ape and ggnewscale, all
  suggested.
* The Marine Honour example data now hold the whole catalogue (281 MAGs) and the
  tree of the 280 bacterial ones (`bac120.tree`).

# magreport 0.5.1

* `plot_abundance_timeline()` gains `condition_labels` and `condition_title`, to
  show readable names in the legend instead of the codes of the sample names
  (for example "H (high)" instead of "O").

# magreport 0.5.0

## New features

* `read_mimag_quality()` reads `mimag_quality.tsv` (mag-pipeline >= 1.8.0): MIMAG
  category, rRNA and tRNA genes, and the GTDB-Tk classification split into
  ranks. The result works with `count_quality()`, `plot_mag_quality()` and
  `plot_taxonomy()`.
* `read_hydrocarbon_genes()` reads `hydrocarbon_per_mag.tsv` (mag-pipeline >=
  1.9.0), the CANT-HYD and HADEG markers of each MAG, and
  `hydrocarbon_groups()` groups the CANT-HYD markers by process.
* `read_mag_abundance()` reads `mag_abundance.tsv` and adds the relative
  abundance of each MAG in each sample. `parse_sample_names()` turns sample
  names into design columns.
* `plot_abundance_timeline()` follows the abundance of groups of MAGs over time,
  per condition, with the mean and standard deviation of the replicates.
* `plot_hydrocarbon_genes()` draws the markers carried by each MAG.
* `read_pipeline()` also finds the MIMAG, hydrocarbon and abundance tables.
* New example data from the Marine Honour oil spill (PRJNA1263040), produced by
  mag-pipeline v1.9.1.

## Bug fixes

* `hydrocarbon_kos()` listed K10944 (`amoA`/`pmoA`, ammonia or methane
  monooxygenase) as an alkane monooxygenase, so `plot_pathway_presence()`
  counted nitrifiers and methanotrophs as alkane degraders. It is removed, as is
  K18089, which was listed as CYP153: KEGG has no specific orthologue for
  CYP153.
