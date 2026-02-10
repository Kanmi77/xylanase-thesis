#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"
IN="$BASE/data/curated/splits_fasta"
OUT="$BASE/results/alignments"
LOG="$BASE/logs/11_run_mafft_split_alignments.log"

mkdir -p "$OUT"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting MAFFT split alignments" | tee -a "$LOG"

for f in fungal_GH10 fungal_GH11 bacterial_GH10 bacterial_GH11; do
  INFA="$IN/${f}.fasta"
  OUTFA="$OUT/${f}.aln.fasta"

  if [[ ! -s "$INFA" ]]; then
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] SKIP empty: $INFA" | tee -a "$LOG"
    continue
  fi

  echo "[`date '+%Y-%m-%d %H:%M:%S'`] Aligning $INFA -> $OUTFA" | tee -a "$LOG"
  mafft --auto --thread 2 "$INFA" > "$OUTFA"
done

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Done." | tee -a "$LOG"
