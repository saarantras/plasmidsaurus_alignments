# Human RNA-cDNA explainability (transcriptome-only)

Date: 2026-03-03

## Setup
- Reads: `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq`
- Reference: GENCODE human transcripts `v48`
- Index: `data/human_refs/v48/human_transcripts.mmi`
- Calling threshold (default/moderate):
  - `best_query_cov >= 0.50`
  - `best_identity >= 0.80`

## Default result
From `human_cdna_explainability.summary.tsv`:
- Total reads: `4586`
- Reads mapped to human transcriptome (any alignment): `2441` (`53.23%`)
- Reads explained as human RNA cDNA: `620` (`13.52%`)
- Total read bases: `2,527,538`
- Bases in explained reads: `358,293` (`14.18%`)

Interpretation:
- About half the reads have at least some transcriptome alignment signal.
- A smaller subset (`~13.5%` of reads; `~14.2%` of bases) passes moderate coverage+identity criteria and is explainable as human RNA-cDNA-like sequence.

## Threshold sensitivity check
Validation runs:
- Strict (`cov>=0.60`, `id>=0.85`): `347/4586` (`7.57%`)
- Permissive (`cov>=0.40`, `id>=0.75`): `1092/4586` (`23.81%`)

As expected, stricter thresholds reduce explained fraction and permissive thresholds increase it.

## Companion outputs
- Per-read table: `human_cdna_explainability.per_read.tsv` (optional detailed artifact)
- Length-binned summary: `human_cdna_explainability.length_binned.tsv`
- Target table: `human_cdna_explainability.top_targets.tsv`
- 2D threshold grid: `human_cdna_explainability.threshold_grid.tsv`
- Distributions:
  - `human_cdna_explainability.identity_hist.tsv`
  - `human_cdna_explainability.query_cov_hist.tsv`
