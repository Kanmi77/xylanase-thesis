#!/usr/bin/env python3
"""
Build a structure manifest for downloaded PDB/mmCIF files.

Input:
  data/raw/pdb/  (*.cif or *.pdb)

Output:
  results/structures/structure_manifest.csv

Rule:
- For each structure, select the longest protein chain as representative.
"""

import os
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder

BASE = os.path.expanduser("~/xylanase-thesis")
PDB_DIR = os.path.join(BASE, "data/raw/pdb")
OUT = os.path.join(BASE, "results/structures/structure_manifest.csv")

def parse_structure(path: str):
    pid = os.path.basename(path).split(".")[0].upper()
    if path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(pid, path)
    model = next(structure.get_models())

    ppb = PPBuilder()
    chain_info = []

    for chain in model:
        # build polypeptides (protein-only)
        pps = ppb.build_peptides(chain)
        length = sum(len(pp.get_sequence()) for pp in pps)
        if length > 0:
            chain_info.append((chain.id, length))

    # choose longest protein chain
    if not chain_info:
        return pid, "", 0

    chain_info.sort(key=lambda x: x[1], reverse=True)
    return pid, chain_info[0][0], chain_info[0][1]

def main() -> int:
    rows = []
    for fn in sorted(os.listdir(PDB_DIR)):
        if not (fn.endswith(".cif") or fn.endswith(".pdb")):
            continue
        path = os.path.join(PDB_DIR, fn)
        pid, chain_id, chain_len = parse_structure(path)
        rows.append({
            "pdb_id": pid,
            "file_path": path,
            "chosen_chain": chain_id,
            "chain_length": chain_len,
            "file_type": "cif" if fn.endswith(".cif") else "pdb",
        })

    df = pd.DataFrame(rows).sort_values(["pdb_id"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(df)}")
    print("Chains selected (non-empty):", (df["chain_length"] > 0).sum())
    print("Structures with no protein chain parsed:", (df["chain_length"] == 0).sum())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
