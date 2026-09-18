# Example output

These files come from a real `mag-pipeline` run on the small public test dataset
(two metagenomic samples, co-assembly), so their format matches the real output.

Two caveats:

- `coassembly.depth.txt` was truncated to the first 120 contigs to keep the
  package small.
- The CheckM2 completeness and contamination values in `mag_quality.tsv`, and
  its taxonomic ranks, were produced by stand-ins for CheckM2 and GTDB-Tk, so
  they are plausible in format but meaningless biologically. They exist so the
  examples and tests run without a 3 GB and a 94 GB reference database.
