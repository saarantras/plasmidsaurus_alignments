# Human-cDNA Unexplained Reads: Consensus Assessment

Date: 2026-03-03

## Question
After calling human RNA-cDNA explained reads (`query_cov >= 0.50` and `identity >= 0.80`), are the remaining unexplained reads mostly one sequence family, or a heterogeneous mixture?

## Inputs
- Source read set: `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq`
- Human-cDNA per-read calls: `results/LFFHG8_1_pcr1_sub5/human_cdna_explainability.per_read.tsv`
- Unexplained subset extracted to:
  - `results/LFFHG8_1_pcr1_sub5/human_cdna_unexplained.read_ids.txt`
  - `results/LFFHG8_1_pcr1_sub5/human_cdna_unexplained.fastq`

Counts:
- Total reads: `4586`
- Human-cDNA explained reads: `620`
- Unexplained reads: `3966`

## Consensus-from-unexplained workflow
1. Built seeded consensus from unexplained-only reads:
   - Consensus FASTA: `results/LFFHG8_1_pcr1_sub5/LFFHG8_1_pcr1_sub5_human_unexplained.consensus.fasta`
   - Consensus summary: `results/LFFHG8_1_pcr1_sub5/LFFHG8_1_pcr1_sub5_human_unexplained.consensus_summary.tsv`
2. Mapped unexplained reads back to that consensus and quantified best-hit coverage/identity:
   - `results/LFFHG8_1_pcr1_sub5/human_cdna_unexplained_vs_consensus.summary.tsv`
   - `results/LFFHG8_1_pcr1_sub5/human_cdna_unexplained_vs_consensus.threshold_grid.tsv`
3. Assessed overlap-graph component structure from all-vs-all PAF:
   - `results/LFFHG8_1_pcr1_sub5/human_cdna_unexplained_overlap_components.summary.tsv`
   - `results/LFFHG8_1_pcr1_sub5/human_cdna_unexplained_overlap_components.component_sizes.tsv`

## Results
### A) How well unexplained reads are captured by one consensus
From `human_cdna_unexplained_vs_consensus.summary.tsv`:
- Reads with any alignment to consensus: `402 / 3966` (`10.14%`)
- Reads passing moderate threshold (`cov>=0.50`, `id>=0.80`): `60 / 3966` (`1.51%`)
- Bases passing moderate threshold: `1.49%`

Threshold sensitivity:
- `cov>=0.30`, `id>=0.70`: `194 / 3966` (`4.89%`)
- `cov>=0.40`, `id>=0.75`: `115 / 3966` (`2.90%`)
- `cov>=0.50`, `id>=0.80`: `60 / 3966` (`1.51%`)

### B) Overlap-graph heterogeneity
Using all-vs-all overlaps filtered at `min_overlap=150` and `min_identity=0.85`:
- Total unexplained reads evaluated: `3966`
- Reads with no overlap record at all: `773`
- Connected components: `2942`
- Singleton components: `2909`
- Largest component: `893` reads (`22.52%` of unexplained set)
- Next component sizes: `75`, `16`, then many tiny clusters

## Interpretation
The unexplained subset is **heterogeneous**, not well modeled by a single consensus sequence.  
At best, one dominant cluster explains roughly one-fifth to one-quarter of unexplained reads, while most reads are singletons or very small components.
