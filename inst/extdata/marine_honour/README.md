# Marine Honour example

Real output of mag-pipeline v1.9.1 on the public metagenomes of the 2024 Marine
Honour oil spill (BioProject PRJNA1263040): 48 intertidal sediment samples, 8
sampling times from 2 to 185 days after the oil arrived, a low-impact (C) and a
high-impact (O) site, 3 replicates each.

- `mimag_quality.tsv`: the 281 dereplicated MAGs, with their MIMAG quality and
  GTDB-Tk (r232) classification.
- `hydrocarbon_per_mag.tsv`: their CANT-HYD and HADEG markers.
- `mag_abundance.tsv.gz`: mean depth of the 281 MAGs in the 48 samples.
- `bac120.tree`: tree of the 280 bacterial MAGs, from the GTDB-Tk alignment
  of 120 marker genes of all of them (`gtdbtk identify`, `gtdbtk align
  --skip_gtdb_refs`) and `gtdbtk infer` (FastTree, WAG + gamma). Unrooted.

Sample names follow `<time>-<site>-<condition>-<replicate>`. The times T0 to T7
are 2, 17, 33, 48, 63, 93, 124 and 185 days after the oil arrived.
