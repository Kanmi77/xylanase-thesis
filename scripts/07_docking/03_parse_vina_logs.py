#!/usr/bin/env python3
"""
Parse Vina logs for top-15 docking runs and output a tidy CSV.

Inputs:
  ~/xylanase-thesis/docking/logs/vina/*.log
  ~/xylanase-thesis/docking/manifests/dock_manifest_top15.csv

Output:
  ~/xylanase-thesis/results/docking/vina_top15_scores.csv
"""

import os
import re
import csv
import glob
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
LOG_DIR = os.path.join(BASE, "docking/logs/vina")
MANIFEST = os.path.join(BASE, "docking/manifests/dock_manifest_top15.csv")
OUT_DIR = os.path.join(BASE, "results/docking")
OUT_CSV = os.path.join(OUT_DIR, "vina_top15_scores.csv")

# Vina affinity table lines look like:
#   1          -7.5      0.0      0.0
ROW_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$")

def parse_one_log(path: str):
    pdb_tag = os.path.basename(path).split("__")[0]
    lig = os.path.basename(path).split("__", 1)[1].replace(".log", "")
    best_aff = None
    best_mode = None
    rmsd_lb = None
    rmsd_ub = None

    in_table = False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "-----+------------+----------+----------" in line:
                in_table = True
                continue
            if in_table:
                m = ROW_RE.match(line)
                if m:
                    mode = int(m.group(1))
                    aff = float(m.group(2))
                    lb = float(m.group(3))
                    ub = float(m.group(4))
                    if best_aff is None or aff < best_aff:
                        best_aff = aff
                        best_mode = mode
                        rmsd_lb = lb
                        rmsd_ub = ub
                # table ends when blank line after values often appears
    return {
        "pdb_tag": pdb_tag,
        "ligand": lig,
        "best_mode": best_mode,
        "best_affinity_kcal_mol": best_aff,
        "rmsd_lb": rmsd_lb,
        "rmsd_ub": rmsd_ub,
        "log_file": os.path.basename(path),
    }

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load manifest metadata
    manifest = pd.read_csv(MANIFEST, dtype=str).fillna("")
    meta = manifest.set_index("pdb_tag").to_dict(orient="index")

    logs = sorted(glob.glob(os.path.join(LOG_DIR, "*.log")))
    if not logs:
        raise SystemExit(f"ERROR: No Vina log files found in {LOG_DIR}")

    rows = []
    for lp in logs:
        r = parse_one_log(lp)
        # attach manifest columns (gh_family, organism_type, etc.)
        m = meta.get(r["pdb_tag"], {})
        r.update({
            "gh_family": m.get("gh_family",""),
            "pdb_id": m.get("pdb_id",""),
            "chosen_chain": m.get("chosen_chain",""),
            "organism_type": m.get("organism_type",""),
            "organism": m.get("organism",""),
            "uniprot_accession": m.get("uniprot_accession",""),
        })
        rows.append(r)

    out = pd.DataFrame(rows)
    # Sort by best affinity (more negative is better)
    out = out.sort_values(["ligand", "best_affinity_kcal_mol"], ascending=[True, True])
    out.to_csv(OUT_CSV, index=False)
    print(f"Wrote: {OUT_CSV} rows={len(out)}")

if __name__ == "__main__":
    main()
