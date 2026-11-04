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
