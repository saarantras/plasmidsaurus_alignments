#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  scripts/fetch_human_transcriptome_refs.sh [--release v48]

Downloads and caches a pinned human transcriptome reference for transcriptome-only
alignment workflows, then builds a minimap2 index.

Outputs:
  data/human_refs/<release>/human_transcripts.fa.gz
  data/human_refs/<release>/human_transcripts.fa
  data/human_refs/<release>/human_transcripts.mmi
  data/human_refs/<release>/REFERENCE_INFO.tsv
USAGE
}

RELEASE="v48"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --release)
      RELEASE="$2"
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

if ! command -v curl >/dev/null 2>&1; then
  echo "Missing required executable: curl" >&2
  exit 1
fi
if ! command -v md5sum >/dev/null 2>&1; then
  echo "Missing required executable: md5sum" >&2
  exit 1
fi
if ! command -v minimap2 >/dev/null 2>&1; then
  echo "Missing required executable: minimap2" >&2
  exit 1
fi
if ! command -v gzip >/dev/null 2>&1; then
  echo "Missing required executable: gzip" >&2
  exit 1
fi

case "$RELEASE" in
  v48)
    BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48"
    SRC_BASENAME="gencode.v48.transcripts.fa.gz"
    EXPECTED_MD5="e4a4d396cca5dd6d0889248b9e93b42a"
    ;;
  *)
    echo "Unsupported release: $RELEASE" >&2
    echo "Currently supported: v48" >&2
    exit 1
    ;;
esac

OUT_DIR="data/human_refs/${RELEASE}"
mkdir -p "$OUT_DIR"

OUT_GZ="${OUT_DIR}/human_transcripts.fa.gz"
OUT_FA="${OUT_DIR}/human_transcripts.fa"
OUT_MMI="${OUT_DIR}/human_transcripts.mmi"
OUT_INFO="${OUT_DIR}/REFERENCE_INFO.tsv"
SRC_URL="${BASE_URL}/${SRC_BASENAME}"

echo "[1/4] Ensure transcript FASTA.GZ is present"
need_download=1
if [[ -f "$OUT_GZ" ]]; then
  existing_md5="$(md5sum "$OUT_GZ" | awk '{print $1}')"
  if [[ "$existing_md5" == "$EXPECTED_MD5" ]]; then
    need_download=0
    echo "  Cache hit (checksum OK): $OUT_GZ"
  else
    echo "  Existing file checksum mismatch; re-downloading."
  fi
fi

if [[ "$need_download" -eq 1 ]]; then
  tmp="${OUT_GZ}.tmp"
  rm -f "$tmp"
  if ! curl -L --fail --retry 3 --retry-delay 2 -o "$tmp" "$SRC_URL"; then
    rm -f "$tmp"
    cat >&2 <<EOF
Failed to download transcriptome reference from:
  $SRC_URL

Remediation:
  1) Check network access to ftp.ebi.ac.uk and retry.
  2) Or place a cached file at:
     $OUT_GZ
     with md5: $EXPECTED_MD5
EOF
    exit 1
  fi
  got_md5="$(md5sum "$tmp" | awk '{print $1}')"
  if [[ "$got_md5" != "$EXPECTED_MD5" ]]; then
    echo "Checksum mismatch for downloaded transcriptome." >&2
    echo "Expected: $EXPECTED_MD5" >&2
    echo "Got:      $got_md5" >&2
    rm -f "$tmp"
    exit 1
  fi
  mv "$tmp" "$OUT_GZ"
fi

echo "[2/4] Ensure uncompressed FASTA is present"
if [[ ! -f "$OUT_FA" ]]; then
  gzip -dc "$OUT_GZ" > "$OUT_FA"
else
  echo "  Cache hit: $OUT_FA"
fi

echo "[3/4] Ensure minimap2 index is present"
if [[ ! -f "$OUT_MMI" ]]; then
  minimap2 -d "$OUT_MMI" "$OUT_FA"
else
  echo "  Cache hit: $OUT_MMI"
fi

echo "[4/4] Write REFERENCE_INFO.tsv"
download_date="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
actual_md5="$(md5sum "$OUT_GZ" | awk '{print $1}')"
{
  echo -e "key\tvalue"
  echo -e "release\t${RELEASE}"
  echo -e "source_url\t${SRC_URL}"
  echo -e "source_basename\t${SRC_BASENAME}"
  echo -e "expected_md5\t${EXPECTED_MD5}"
  echo -e "actual_md5\t${actual_md5}"
  echo -e "download_or_verify_utc\t${download_date}"
  echo -e "fasta_gz\t${OUT_GZ}"
  echo -e "fasta\t${OUT_FA}"
  echo -e "minimap2_index\t${OUT_MMI}"
} > "$OUT_INFO"

echo "Done."
echo "Transcript FASTA:  $OUT_FA"
echo "Transcript index:  $OUT_MMI"
echo "Reference info:    $OUT_INFO"
