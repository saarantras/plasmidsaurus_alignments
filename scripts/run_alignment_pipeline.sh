#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  scripts/run_alignment_pipeline.sh \
    --sample <sample_id> \
    --reads <reads.fastq.gz> \
    --reference-gb <reference.gb|reference.gbk> \
    [--ref-name <reference_name>] \
    [--threads <n>]

Outputs are written under results/<sample_id>/.
USAGE
}

SAMPLE_ID=""
READS_FASTQ=""
REFERENCE_GB=""
REF_NAME=""
THREADS="4"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)
      SAMPLE_ID="$2"
      shift 2
      ;;
    --reads)
      READS_FASTQ="$2"
      shift 2
      ;;
    --reference-gb)
      REFERENCE_GB="$2"
      shift 2
      ;;
    --ref-name)
      REF_NAME="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$SAMPLE_ID" || -z "$READS_FASTQ" || -z "$REFERENCE_GB" ]]; then
  usage
  exit 1
fi

for tool in minimap2 samtools bedtools python; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Missing required executable in PATH: $tool" >&2
    exit 1
  fi
done

if [[ ! -f "$READS_FASTQ" ]]; then
  echo "Reads file not found: $READS_FASTQ" >&2
  exit 1
fi

if [[ ! -f "$REFERENCE_GB" ]]; then
  echo "Reference GenBank not found: $REFERENCE_GB" >&2
  exit 1
fi

if [[ -z "$REF_NAME" ]]; then
  basename_ref="$(basename "$REFERENCE_GB")"
  REF_NAME="${basename_ref%.*}"
fi

REF_FASTA="data/references/${REF_NAME}.fasta"
REF_BED="data/references/${REF_NAME}.features.bed"
REF_INDEX="data/references/${REF_NAME}.mmi"
OUT_DIR="results/${SAMPLE_ID}"
BAM_PATH="${OUT_DIR}/${SAMPLE_ID}.bam"
DEPTH_PATH="${OUT_DIR}/${SAMPLE_ID}.depth.tsv"
FEATURE_COVERAGE_RAW="${OUT_DIR}/${SAMPLE_ID}.feature_coverage.raw.tsv"
FEATURE_COVERAGE_SUMMARY="${OUT_DIR}/${SAMPLE_ID}.feature_coverage.summary.tsv"
DEPTH_PNG="${OUT_DIR}/${SAMPLE_ID}.depth_with_features.png"
FEATURE_PNG="${OUT_DIR}/${SAMPLE_ID}.feature_mean_depth.png"

mkdir -p "$(dirname "$REF_FASTA")" "$OUT_DIR"

echo "[1/7] Convert GenBank to FASTA + BED"
python scripts/gb_to_fasta_and_bed.py \
  --genbank "$REFERENCE_GB" \
  --fasta "$REF_FASTA" \
  --bed "$REF_BED"

echo "[2/7] Build minimap2 index"
minimap2 -d "$REF_INDEX" "$REF_FASTA"

echo "[3/7] Align reads -> sorted BAM"
minimap2 -ax map-ont -t "$THREADS" "$REF_INDEX" "$READS_FASTQ" \
  | samtools sort -@ "$THREADS" -o "$BAM_PATH"

echo "[4/7] Index BAM"
samtools index "$BAM_PATH"

echo "[5/7] Per-base depth"
samtools depth -aa "$BAM_PATH" > "$DEPTH_PATH"

echo "[6/7] Feature coverage with bedtools"
bedtools coverage \
  -a "$REF_BED" \
  -b "$BAM_PATH" \
  > "$FEATURE_COVERAGE_RAW"

echo "[7/7] Feature summary + plots"
python scripts/summarize_feature_coverage.py \
  --features "$REF_BED" \
  --depth "$DEPTH_PATH" \
  --output "$FEATURE_COVERAGE_SUMMARY" \
  --plot "$FEATURE_PNG"

python scripts/plot_depth_with_features.py \
  --depth "$DEPTH_PATH" \
  --features "$REF_BED" \
  --output "$DEPTH_PNG" \
  --title "${SAMPLE_ID}: depth with features"

echo "Done. Outputs in $OUT_DIR"
