# Marine Honour example

Real output of mag-pipeline v1.9.1 on the public metagenomes of the 2024 Marine
Honour oil spill (BioProject PRJNA1263040): 48 intertidal sediment samples, 8
sampling times from 2 to 185 days after the oil arrived, a low-impact (C) and a
high-impact (O) site, 3 replicates each.

- `mimag_quality.tsv`: 74 of the 281 dereplicated MAGs (the 62 with a CANT-HYD
  marker, the 6 of high MIMAG quality and the strongest responders to the
  spill).
- `hydrocarbon_per_mag.tsv`: CANT-HYD and HADEG markers of those 74 MAGs.
- `mag_abundance.tsv.gz`: mean depth of all 281 MAGs in the 48 samples, so that
  relative abundances are computed over the whole catalogue.

Sample names follow `<time>-<site>-<condition>-<replicate>`. The times T0 to T7
are 2, 17, 33, 48, 63, 93, 124 and 185 days after the oil arrived.
