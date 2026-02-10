#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"
INPUTS="$BASE/data/curated/foldx_inputs.csv"
PDBDIR="$BASE/foldx/pdb_chains"
OUTDIR="$BASE/foldx/repair"
LOG="$BASE/logs/18B_foldx_repairpdb.log"

# ---- SET THIS ----
FOLDX_BIN="/home/ubuntu/xylanase-thesis/foldx/bin/FoldX"
# If your binary is different, change it above.

mkdir -p "$OUTDIR"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] FoldX RepairPDB start" | tee -a "$LOG"
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Using: $FOLDX_BIN" | tee -a "$LOG"

# sanity
if [[ ! -x "$FOLDX_BIN" ]]; then
  echo "ERROR: FoldX binary not executable at: $FOLDX_BIN" | tee -a "$LOG"
  exit 1
fi

# Run RepairPDB per file
tail -n +2 "$INPUTS" | cut -d',' -f1,3,0 >/dev/null 2>&1 || true

# We will read the foldx_pdb column by header name using python to avoid CSV parsing issues
python - <<'PY'
import pandas as pd, os, subprocess, shlex
BASE=os.path.expanduser("~/xylanase-thesis")
inputs=os.path.join(BASE,"data/curated/foldx_inputs.csv")
outdir=os.path.join(BASE,"foldx/repair")
foldx=os.path.join(BASE,"foldx/foldx_5/foldx5_Linux64/FoldX")

df=pd.read_csv(inputs, dtype=str).fillna("")
pdbs=df["foldx_pdb"].tolist()

for p in pdbs:
    if not os.path.exists(p): 
        continue
    # FoldX writes outputs to current working dir, so run inside outdir
    cmd=[foldx, "--command=RepairPDB", f"--pdb={os.path.basename(p)}"]
    # copy input PDB into outdir for foldx to find it
    dst=os.path.join(outdir, os.path.basename(p))
    if not os.path.exists(dst):
        import shutil; shutil.copy2(p,dst)
    subprocess.run(cmd, cwd=outdir, check=False)
print("DONE")
PY

echo "[`date '+%Y-%m-%d %H:%M:%S'`] FoldX RepairPDB done" | tee -a "$LOG"
