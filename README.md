# magreport

<!-- badges: start -->
[![R-CMD-check](https://github.com/AntonioPuriel/magreport/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AntonioPuriel/magreport/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

R package to explore and report MAG recovery results from
[mag-pipeline](https://github.com/AntonioPuriel/mag-pipeline).

The pipeline writes plain tables; `magreport` reads them into tidy data frames
with the right types, so you can go straight to analysis and figures.

> 🚧 Work in progress — readers and plots are implemented; the HTML report is next.

## Installation

```r
# install.packages("remotes")
remotes::install_github("AntonioPuriel/magreport")
```

## Usage

```r
library(magreport)

results <- "path/to/mag-pipeline/results"

mags    <- read_mag_quality(file.path(results, "05_quality/mag_quality.tsv"))
bins    <- read_bin_summary(file.path(results, "04_binning/bin_summary.tsv"))
mapping <- read_mapping_summary(file.path(results, "03_mapping/mapping_summary.tsv"))
depth   <- read_depth(file.path(results, "03_mapping/coassembly.depth.txt"))

count_quality(mags, by = "binner")

plot_mag_quality(mags, facet_by = "binner")
plot_mapping_rates(mapping)
plot_taxonomy(mags, rank = "phylum", fill_by = "quality")
```

The plotting functions need ggplot2, which is a suggested dependency:
`install.packages("ggplot2")`.

Every function ships an example that runs on the output included in the package,
so you can try it without your own data:

```r
path <- system.file("extdata", "example_run", "mag_quality.tsv", package = "magreport")
head(read_mag_quality(path))
```

## Functions

| Function | Reads |
|----------|-------|
| `read_mag_quality()` | `05_quality/mag_quality.tsv` (quality and taxonomy per MAG) |
| `read_bin_summary()` | `04_binning/bin_summary.tsv` (size, N50 and GC per bin) |
| `read_mapping_summary()` | `03_mapping/mapping_summary.tsv` (alignment rate per sample) |
| `read_depth()` | `03_mapping/*.depth.txt` (contig coverage per sample) |
| `count_quality()` | Counts MAGs per MIMAG quality category |

## Plots

All figures below come from the example output shipped with the package, and are
regenerated with `data-raw/make_readme_figures.R`.

### `plot_mag_quality()`

Completeness against contamination, with the MIMAG thresholds and the
high-quality corner shaded. Point size is the size of the MAG.

![MAG quality](man/figures/README-mag-quality.png)

### `plot_mapping_rates()`

Share of reads of each sample that map back to the assembly. Samples below the
threshold are highlighted, since their coverage is less trustworthy.

![Mapping rates](man/figures/README-mapping-rates.png)

### `plot_taxonomy()`

MAGs per taxon at any rank. Unclassified MAGs are kept as their own group
instead of being dropped.

![Taxonomy](man/figures/README-taxonomy.png)

## Documentation

Full function reference: <https://antoniopuriel.github.io/magreport/>

## License

MIT
