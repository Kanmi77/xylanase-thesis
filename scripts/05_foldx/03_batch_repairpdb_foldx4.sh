#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"
FOLDX="$BASE/foldx/bin/FoldX"
INPUT_DIR="$BASE/foldx/pdb_chains"
OUTDIR="$BASE/foldx/repair"
LOG="$BASE/logs/18B_batch_repairpdb_foldx4.log"
ROTABASE_SRC="$BASE/foldx/bin/rotabase.txt"

mkdir -p "$OUTDIR"
mkdir -p "$(dirname "$LOG")"

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Batch RepairPDB start" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] FoldX: $FOLDX" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Input dir: $INPUT_DIR" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Output dir: $OUTDIR" | tee -a "$LOG"

[[ -x "$FOLDX" ]] || { echo "ERROR: FoldX not executable: $FOLDX" | tee -a "$LOG"; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: Missing input dir: $INPUT_DIR" | tee -a "$LOG"; exit 1; }
[[ -s "$ROTABASE_SRC" ]] || { echo "ERROR: Missing rotabase: $ROTABASE_SRC" | tee -a "$LOG"; exit 1; }

# Put rotabase in the *working directory* used for FoldX
cp -f "$ROTABASE_SRC" "$OUTDIR/rotabase.txt"

TOTAL=$(ls -1 "$INPUT_DIR"/*.pdb 2>/dev/null | wc -l)
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Found input PDBs: $TOTAL" | tee -a "$LOG"

OK=0
SKIP=0
FAIL=0
FAIL_RC=0
FAIL_NOFILE=0

for P in "$INPUT_DIR"/*.pdb; do
  BN=$(basename "$P")
  TAG="${BN%.pdb}"
  REPAIRED="$OUTDIR/${TAG}_Repair.pdb"

  if [[ -s "$REPAIRED" ]]; then
    SKIP=$((SKIP+1))
    continue
  fi

  # Copy the input PDB into OUTDIR (FoldX must see it in its CWD)
  cp -f "$P" "$OUTDIR/$BN"

  echo "[`date '+%Y-%m-%d %H:%M:%S'`] Repairing $BN" | tee -a "$LOG"

  set +e
  ( cd "$OUTDIR" && "$FOLDX" --command=RepairPDB --pdb="$BN" ) > "$OUTDIR/${BN}.repair.stdout.txt" 2>&1
  RC=$?
  set -e

  if [[ $RC -ne 0 ]]; then
    FAIL=$((FAIL+1))
    FAIL_RC=$((FAIL_RC+1))
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] FAIL rc=$RC $BN" | tee -a "$LOG"
    continue
  fi

  if [[ -s "$REPAIRED" ]]; then
    OK=$((OK+1))
  else
    FAIL=$((FAIL+1))
    FAIL_NOFILE=$((FAIL_NOFILE+1))
    echo "[`date '+%Y-%m-%d %H:%M:%S'`] FAIL no repaired output for $BN (expected $REPAIRED)" | tee -a "$LOG"
  fi
done

echo "[`date '+%Y-%m-%d %H:%M:%S'`] Done. ok=$OK skip=$SKIP fail=$FAIL (rc_fail=$FAIL_RC nofile_fail=$FAIL_NOFILE)" | tee -a "$LOG"
