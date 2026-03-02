# LFFHG8 vs gblock_f1r1 conclusions (March 2, 2026)

## Summary
- The expected amplicon product is present, but at very low abundance relative to other sequence content.
- Most reads are non-target/background relative to `gblock_f1r1`.

## Evidence
1. Whole dataset mapping to reference (`gblock_f1r1`):
   - `4 / 4586` reads mapped (`0.09%`).
   - Source: `results/LFFHG8_1_pcr1_sub5/LFFHG8_1_pcr1_sub5.flagstat.txt`.

2. Forward-primer check (`F1 primer = AATGATAATAATCGCGTCGACG`):
   - Exactly `3` reads contain the forward primer sequence.
   - In all 3 reads, the primer starts at read position 0 (near the read start as expected).

3. Alignment quality for the 3 forward-primer-positive reads:
   - `795_LFFHG8_1`: 98.64% identity, nearly full-length match to reference.
   - `3446_LFFHG8_1`: 97.00% identity, long match to reference.
   - `3685_LFFHG8_1`: partial/lower-quality match.
   - All three map to `gblock_f1r1` with MAPQ 60.

4. Consensus from dominant non-reference cluster:
   - Reference-free seeded consensus length: 1611 bp.
   - Built from a seed-cluster with `1046` mapped reads (21.54% of reads to seed).
   - Consensus vs `gblock_f1r1` overlap is limited to a 27 bp exact region at the reference end.

5. Feature overlap from consensus vs reference:
   - The 27 bp exact overlap corresponds only to end-primer region features (`R1 primer` / `R3 primer`).
   - No broad overlap across the full annotated construct.

## Interpretation
- Data strongly suggest true target/product exists in a small read subset.
- The majority of sequencing content is off-target/background ("gorp").
- Off-target sequence identification is in progress via BLAST of the non-reference consensus.
