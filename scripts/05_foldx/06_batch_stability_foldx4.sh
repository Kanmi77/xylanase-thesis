#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"
FOLDX="$BASE/foldx/bin/FoldX"
REPAIR_DIR="$BASE/foldx/repair"
OUTDIR="$BASE/foldx/stability"
LOG="$BASE/logs/19A_batch_stability_foldx4.log"

mkdir -p "$OUTDIR"
mkdir -p "$(dirname "$LOG")"

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Batch Stability start" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] FoldX: $FOLDX" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Repair dir: $REPAIR_DIR" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Output dir: $OUTDIR" | tee -a "$LOG"

[[ -x "$FOLDX" ]] || { echo "ERROR: FoldX not executable: $FOLDX" | tee -a "$LOG"; exit 1; }
[[ -d "$REPAIR_DIR" ]] || { echo "ERROR: Missing repair dir: $REPAIR_DIR" | tee -a "$LOG"; exit 1; }

# Ensure rotabase exists in OUTDIR working directory for FoldX
if [[ -s "$REPAIR_DIR/rotabase.txt" ]]; then
  cp -f "$REPAIR_DIR/rotabase.txt" "$OUTDIR/rotabase.txt"
elif [[ -s "$BASE/foldx/bin/rotabase.txt" ]]; then
  cp -f "$BASE/foldx/bin/rotabase.txt" "$OUTDIR/rotabase.txt"
else
  echo "ERROR: rotabase.txt not found in repair or bin" | tee -a "$LOG"
  exit 1
fi

TOTAL=$(ls -1 "$REPAIR_DIR"/*_Repair.pdb 2>/dev/null | wc -l)
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Found repaired PDBs: $TOTAL" | tee -a "$LOG"

OK=0
SKIP=0
FAIL=0

for P in "$REPAIR_DIR"/*_Repair.pdb; do
  BN=$(basename "$P")

  # FoldX usually produces fxout; we keep a per-structure log too
  OUT_STD="$OUTDIR/${BN}.stability.stdout.txt"

  # If we already ran Stability for this BN, skip
  if [[ -s "$OUT_STD" ]]; then
    SKIP=$((SKIP+1))
    continue
  fi

  # Copy repaired PDB into OUTDIR working directory
  cp -f "$P" "$OUTDIR/$BN"

  echo "[`date '+%Y-%m-%d %H:%M:%S'`] Stability $BN" | tee -a "$LOG"
  set +e
  ( cd "$OUTDIR" && "$FOLDX" --command=Stability --pdb="$BN" ) > "$OUT_STD" 2>&1
  RC=$?
  set -e

  if [[ $RC -ne 0 ]]; then
    FAIL=$((FAIL+1))
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] FAIL rc=$RC $BN" | tee -a "$LOG"
    continue
  fi

  OK=$((OK+1))
done

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Done. ok=$OK skip=$SKIP fail=$FAIL" | tee -a "$LOG"
