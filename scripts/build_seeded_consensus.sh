#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  scripts/build_seeded_consensus.sh \
    --sample <sample_id> \
    --reads <reads.fastq|reads.fastq.gz> \
    [--outdir results/<sample_id>] \
    [--threads 4] \
    [--min-overlap 150] \
    [--min-ident 0.85] \
    [--min-depth 3]

Builds a reference-free consensus by:
1) all-vs-all read overlaps,
2) picking the most connected seed read,
3) mapping reads to seed,
4) calling consensus with samtools consensus.
USAGE
}

SAMPLE_ID=""
READS=""
OUTDIR=""
THREADS="4"
MIN_OVERLAP="150"
MIN_IDENT="0.85"
MIN_DEPTH="3"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)
      SAMPLE_ID="$2"; shift 2 ;;
    --reads)
      READS="$2"; shift 2 ;;
    --outdir)
      OUTDIR="$2"; shift 2 ;;
    --threads)
      THREADS="$2"; shift 2 ;;
    --min-overlap)
      MIN_OVERLAP="$2"; shift 2 ;;
    --min-ident)
      MIN_IDENT="$2"; shift 2 ;;
    --min-depth)
      MIN_DEPTH="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 1 ;;
  esac
done

if [[ -z "$SAMPLE_ID" || -z "$READS" ]]; then
  usage
  exit 1
fi

if [[ -z "$OUTDIR" ]]; then
  OUTDIR="results/${SAMPLE_ID}"
fi

for t in minimap2 samtools awk; do
  if ! command -v "$t" >/dev/null 2>&1; then
    echo "Missing tool in PATH: $t" >&2
    exit 1
  fi
done

if [[ ! -f "$READS" ]]; then
  echo "Reads not found: $READS" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

READS_FASTA="$OUTDIR/${SAMPLE_ID}.reads.fasta"
PAF="$OUTDIR/${SAMPLE_ID}.all_vs_all.paf"
SEED_ID_TXT="$OUTDIR/${SAMPLE_ID}.consensus_seed_id.txt"
SEED_FASTA="$OUTDIR/${SAMPLE_ID}.consensus_seed.fasta"
SEED_BAM="$OUTDIR/${SAMPLE_ID}.consensus_seed_map.bam"
SEED_FLAGSTAT="$OUTDIR/${SAMPLE_ID}.consensus_seed_map.flagstat.txt"
SEED_DEPTH="$OUTDIR/${SAMPLE_ID}.consensus_seed_map.depth.tsv"
CONSENSUS_FASTA="$OUTDIR/${SAMPLE_ID}.consensus.fasta"
SUMMARY_TSV="$OUTDIR/${SAMPLE_ID}.consensus_summary.tsv"

echo "[1/7] Convert reads to FASTA"
if [[ "$READS" == *.gz ]]; then
  gzip -dc "$READS" | awk 'NR%4==1{h=substr($0,2)} NR%4==2{print ">"h"\n"$0}' > "$READS_FASTA"
else
  awk 'NR%4==1{h=substr($0,2)} NR%4==2{print ">"h"\n"$0}' "$READS" > "$READS_FASTA"
fi

echo "[2/7] All-vs-all overlaps"
minimap2 -x ava-ont -t "$THREADS" "$READS_FASTA" "$READS_FASTA" > "$PAF"

echo "[3/7] Pick most connected seed read"
seed_id="$({ awk -v min_ov="$MIN_OVERLAP" -v min_id="$MIN_IDENT" 'BEGIN{FS="\t"} $1!=$6{block=$11+0; id=($10+0)/(($11+0)==0?1:($11+0)); if(block>=min_ov && id>=min_id){c[$1]++}} END{best=""; max=-1; for(k in c){if(c[k]>max){max=c[k]; best=k}}; print best}' "$PAF"; } )"
if [[ -z "$seed_id" ]]; then
  echo "Could not identify seed read. Try lowering --min-overlap or --min-ident." >&2
  exit 1
fi

echo "$seed_id" > "$SEED_ID_TXT"
if [[ "$READS" == *.gz ]]; then
  gzip -dc "$READS" | awk -v target="$seed_id" 'NR%4==1{h=substr($0,2)} NR%4==2{if(h==target){print ">"h"\n"$0; exit}}' > "$SEED_FASTA"
else
  awk -v target="$seed_id" 'NR%4==1{h=substr($0,2)} NR%4==2{if(h==target){print ">"h"\n"$0; exit}}' "$READS" > "$SEED_FASTA"
fi

seed_len=$(awk 'BEGIN{l=0} /^>/{next} {l+=length($0)} END{print l}' "$SEED_FASTA")

echo "[4/7] Map reads to seed"
minimap2 -ax map-ont -t "$THREADS" "$SEED_FASTA" "$READS" \
  | samtools sort -@ "$THREADS" -o "$SEED_BAM"
samtools index "$SEED_BAM"
samtools flagstat "$SEED_BAM" > "$SEED_FLAGSTAT"

echo "[5/7] Depth on seed mapping"
samtools depth -aa "$SEED_BAM" > "$SEED_DEPTH"

echo "[6/7] Call consensus"
samtools consensus -f FASTA -a -d "$MIN_DEPTH" -T "$SEED_FASTA" -o "$CONSENSUS_FASTA" "$SEED_BAM"

cons_len=$(awk 'BEGIN{l=0} /^>/{next} {l+=length($0)} END{print l}' "$CONSENSUS_FASTA")
cons_n=$(awk 'BEGIN{n=0} /^>/{next} {t=toupper($0); gsub(/[^N]/,"",t); n+=length(t)} END{print n}' "$CONSENSUS_FASTA")
mapped_total=$(awk '$4=="mapped"{print $1; exit}' "$SEED_FLAGSTAT")
mapped_pct=$(awk '$4=="mapped"{v=$5; gsub(/[()%]/,"",v); print v; exit}' "$SEED_FLAGSTAT")
map_line="${mapped_total} (${mapped_pct})"

depth_stats=$(awk '{sum+=$3; if($3>0)c++; n++} END{printf "%d\t%d\t%.3f\t%.4f", n, c, sum/n, c/n}' "$SEED_DEPTH")
positions=$(echo "$depth_stats" | cut -f1)
covered=$(echo "$depth_stats" | cut -f2)
mean_depth=$(echo "$depth_stats" | cut -f3)
covered_frac=$(echo "$depth_stats" | cut -f4)

echo "[7/7] Write summary"
{
  echo -e "sample_id\tseed_read_id\tseed_length\tconsensus_length\tconsensus_N_count\tmapped_reads\tpositions\tcovered_positions\tmean_depth\tcovered_fraction"
  echo -e "${SAMPLE_ID}\t${seed_id}\t${seed_len}\t${cons_len}\t${cons_n}\t${map_line}\t${positions}\t${covered}\t${mean_depth}\t${covered_frac}"
} > "$SUMMARY_TSV"

echo "Consensus FASTA: $CONSENSUS_FASTA"
echo "Summary TSV: $SUMMARY_TSV"
