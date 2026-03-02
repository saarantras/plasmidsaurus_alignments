# plasmidsaurus_alignments

Align long reads from Plasmidsaurus to reference sequences (plasmid or amplicon) and generate read-density maps, including coverage over annotated GenBank features.

## Environment
Conda env used: `plasmidsaurus-align`

Environment exports are committed under `env/`:
- `env/environment.plasmidsaurus-align.yml` (portable conda spec)
- `env/environment.plasmidsaurus-align.from-history.yml` (minimal requested specs)
- `env/environment.plasmidsaurus-align.lock.txt` (`conda --explicit` lock-style spec)
- `env/environment.plasmidsaurus-align.pip.txt` (`pip freeze`)

Refresh these files after package changes:
```bash
scripts/export_env.sh plasmidsaurus-align
```

## Implemented pipeline
Driver script:
```bash
conda activate plasmidsaurus-align
scripts/run_alignment_pipeline.sh \
  --sample <sample_id> \
  --reads data/reads/<sample>.fastq.gz \
  --reference-gb data/references/<ref>.gb \
  --threads 4
```

Non-activated shell alternative:
```bash
conda run -n plasmidsaurus-align bash scripts/run_alignment_pipeline.sh \
  --sample <sample_id> \
  --reads data/reads/<sample>.fastq.gz \
  --reference-gb data/references/<ref>.gb \
  --threads 4
```

What it runs:
1. `scripts/gb_to_fasta_and_bed.py` to convert GenBank to FASTA + BED feature intervals.
2. `minimap2` index + alignment to BAM.
3. `samtools depth -aa` for base-wise coverage.
4. `bedtools coverage` for raw per-feature overlap metrics.
5. `scripts/summarize_feature_coverage.py` for feature-level summary table + mean-depth plot.
6. `scripts/plot_depth_with_features.py` for depth line plot with annotation tracks.

Outputs:
- `results/<sample>/<sample>.bam` and `.bam.bai`
- `results/<sample>/<sample>.flagstat.txt`
- `results/<sample>/<sample>.depth.tsv`
- `results/<sample>/<sample>.feature_coverage.raw.tsv`
- `results/<sample>/<sample>.feature_coverage.summary.tsv`
- `results/<sample>/<sample>.depth_with_features.png`
- `results/<sample>/<sample>.feature_mean_depth.png`

## LFFHG8 + gblock_f1r1 run status
- Reference provided and used: `data/references/gblock_f1r1.gb`
- FASTQ used: `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq` (unpacked from `data/reads/LFFHG8_fastq.zip`)
- Alignment run completed successfully.
- `samtools flagstat`: `4 / 4586` reads mapped (`0.09%`).
- Feature coverage and plots were generated in `results/LFFHG8_1_pcr1_sub5/`.
- Earlier mock smoke-test outputs are preserved in `results/LFFHG8_1_pcr1_sub5_mock/`.

## Suggested layout
```text
data/
  reads/
  references/
scripts/
results/
```

## Sample mapping table
Use this table to track which read file is mapped to which reference and whether the target is an amplicon or plasmid.

| sample_id | reads_file | reference_genbank | reference_seq_name | status (amplicon/plasmid) | notes |
|---|---|---|---|---|---|
| LFFHG8_1_pcr1_sub5 | `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq` | `data/references/gblock_f1r1.gb` | `gblock_f1r1` | amplicon | Aligned; 4/4586 reads mapped (0.09%) |
| LFFHG8_1_pcr1_sub5_mock | `data/reads/LFFHG8_1_pcr1_sub5.mock.fastq` | `data/references/gblock_f1r1.gb` | `gblock_f1r1` | amplicon | Historical smoke-test input generated from consensus FASTA |
