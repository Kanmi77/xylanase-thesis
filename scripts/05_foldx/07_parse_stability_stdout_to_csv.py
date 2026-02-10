#!/usr/bin/env python3
"""
Parse FoldX4 Stability stdout logs into a clean CSV.

Input:
  foldx/stability/*_Repair.pdb.stability.stdout.txt

Output:
  results/foldx/foldx4_wt_stability.csv

Extracts:
- pdb_tag (e.g., 1AXK_B_Repair)
- total_energy (float) from 'Total = ...' line
"""

import os, re, glob
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_DIR = os.path.join(BASE, "foldx/stability")
OUT = os.path.join(BASE, "results/foldx/foldx4_wt_stability.csv")

TOTAL_RE = re.compile(r"^\s*Total\s*=\s*([-+]?\d+(?:\.\d+)?)", re.IGNORECASE)

def parse_one(path: str):
    bn = os.path.basename(path)
    # bn looks like: 1AXK_B_Repair.pdb.stability.stdout.txt
    pdb_tag = bn.replace(".pdb.stability.stdout.txt", "")
    energy = None
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = TOTAL_RE.search(line.strip())
            if m:
                try:
                    energy = float(m.group(1))
                except:
                    energy = None
                break
    return pdb_tag, energy

def main() -> int:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    logs = sorted(glob.glob(os.path.join(IN_DIR, "*_Repair.pdb.stability.stdout.txt")))
    if not logs:
        raise SystemExit(f"ERROR: No stability stdout logs found in {IN_DIR}")

    rows = []
    for p in logs:
        tag, e = parse_one(p)
        rows.append({"pdb_tag": tag, "foldx_wt_total_energy": e, "stdout_file": os.path.basename(p)})

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(df)}")
    print("Non-null energies:", df["foldx_wt_total_energy"].notna().sum())
    print("Null energies:", df["foldx_wt_total_energy"].isna().sum())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
