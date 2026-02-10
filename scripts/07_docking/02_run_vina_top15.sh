
#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/xylanase-thesis"

MANIFEST="$BASE/docking/manifests/dock_manifest_top15.csv"
LIG_DIR="$BASE/docking/ligands"
OUT_DIR="$BASE/docking/results"
LOG_DIR="$BASE/docking/logs/vina"
RES_DIR="$BASE/results/docking"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$RES_DIR"

LIG1="$LIG_DIR/xylobiose.pdbqt"
LIG2="$LIG_DIR/xylotriose.pdbqt"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: Missing manifest: $MANIFEST" >&2
  exit 1
fi

if [[ ! -f "$LIG1" || ! -f "$LIG2" ]]; then
  echo "ERROR: Missing ligand pdbqt(s). Expected:" >&2
  echo "  $LIG1" >&2
  echo "  $LIG2" >&2
  exit 1
fi

if ! command -v vina >/dev/null 2>&1; then
  echo "ERROR: vina not in PATH. Activate the env where vina is installed." >&2
  exit 1
fi

echo "[INFO] Manifest: $MANIFEST"
echo "[INFO] Ligands : $LIG1 , $LIG2"
echo "[INFO] Out dir : $OUT_DIR"
echo "[INFO] Log dir : $LOG_DIR"

TMP_TSV="$(mktemp)"
trap 'rm -f "$TMP_TSV"' EXIT

# Emit safe TSV: pdb_tag \t receptor_pdbqt \t vina_conf
python - <<'PY' > "$TMP_TSV"
import csv, os
base=os.path.expanduser("~/xylanase-thesis")
manifest=os.path.join(base,"docking/manifests/dock_manifest_top15.csv")
with open(manifest, newline="") as f:
    r=csv.DictReader(f)
    for row in r:
        print("\t".join([row["pdb_tag"], row["receptor_pdbqt"], row["vina_conf"]]))
PY

# Loop over TSV lines
while IFS=$'\t' read -r pdb_tag receptor_pdbqt vina_conf; do
  if [[ -z "$pdb_tag" ]]; then
    continue
  fi

  if [[ ! -f "$receptor_pdbqt" ]]; then
    echo "WARN: missing receptor pdbqt for $pdb_tag: $receptor_pdbqt (skip)"
    continue
  fi

  if [[ ! -f "$vina_conf" ]]; then
    echo "WARN: missing vina conf for $pdb_tag: $vina_conf (skip)"
    continue
  fi

  for lig in "$LIG1" "$LIG2"; do
    lig_base="$(basename "$lig" .pdbqt)"
    out_pose="$OUT_DIR/${pdb_tag}__${lig_base}.out.pdbqt"
    out_log="$LOG_DIR/${pdb_tag}__${lig_base}.log"

    echo "[RUN] $pdb_tag + $lig_base"
    vina --config "$vina_conf" \
         --ligand "$lig" \
         --out "$out_pose" \
         --log "$out_log"
  done
done < "$TMP_TSV"

echo "[DONE] Vina docking finished."
x
