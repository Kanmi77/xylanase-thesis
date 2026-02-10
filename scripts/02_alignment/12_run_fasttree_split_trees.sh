#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"
IN="$BASE/results/alignments"
OUT="$BASE/results/trees"
LOG="$BASE/logs/12_run_fasttree_split_trees.log"

mkdir -p "$OUT"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Starting FastTree on split alignments" | tee -a "$LOG"

for f in fungal_GH10 fungal_GH11 bacterial_GH10 bacterial_GH11; do
  ALN="$IN/${f}.aln.fasta"
  TREE="$OUT/${f}.fasttree.nwk"

  if [[ ! -s "$ALN" ]]; then
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] SKIP missing/empty: $ALN" | tee -a "$LOG"
    continue
  fi

  echo "[`date '+%Y-%m-%d %H:%M:%S'`] Tree: $ALN -> $TREE" | tee -a "$LOG"
  # protein mode: -lg (LG model) is a solid default; -gamma improves rate heterogeneity
  FastTree -lg -gamma "$ALN" > "$TREE"
done

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Done." | tee -a "$LOG"
