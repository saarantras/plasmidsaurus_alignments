#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-plasmidsaurus-align}"
OUT_DIR="env"

mkdir -p "$OUT_DIR"

YML_PATH="$OUT_DIR/environment.${ENV_NAME}.yml"
HISTORY_YML_PATH="$OUT_DIR/environment.${ENV_NAME}.from-history.yml"
LOCK_PATH="$OUT_DIR/environment.${ENV_NAME}.lock.txt"
PIP_PATH="$OUT_DIR/environment.${ENV_NAME}.pip.txt"

echo "Exporting conda environment: $ENV_NAME"
conda env export -n "$ENV_NAME" --no-builds > "$YML_PATH"
conda env export -n "$ENV_NAME" --from-history > "$HISTORY_YML_PATH"
conda list -n "$ENV_NAME" --explicit > "$LOCK_PATH"
conda run -n "$ENV_NAME" python -m pip freeze > "$PIP_PATH"

echo "Wrote: $YML_PATH"
echo "Wrote: $HISTORY_YML_PATH"
echo "Wrote: $LOCK_PATH"
echo "Wrote: $PIP_PATH"
