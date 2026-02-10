#!/usr/bin/env python3
"""
Build a normalized PDB inventory from the curated UniProt master CSV.

Input:
  data/curated/xylanase_master_uniprot.csv

Output:
  data/curated/pdb_inventory.csv
    columns:
      uniprot_accession, organism_type, organism, pdb_id
"""

import os, re
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")
OUT_CSV = os.path.join(BASE, "data/curated/pdb_inventory.csv")

PDB_RE = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")  # e.g. 3U7B, 1XYZ

def main() -> int:
    df = pd.read_csv(IN_CSV, dtype=str).fillna("")
    rows = []
    for _, r in df.iterrows():
        acc = r.get("uniprot_accession", "").strip()
        pdb_field = r.get("xref_pdb", "")
        pdb_ids = PDB_RE.findall(str(pdb_field))
        for pid in pdb_ids:
            rows.append({
                "uniprot_accession": acc,
                "organism_type": r.get("organism_type", ""),
                "organism": r.get("organism", ""),
                "pdb_id": pid.upper()
            })

    inv = pd.DataFrame(rows).drop_duplicates()
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    inv.to_csv(OUT_CSV, index=False)

    print(f"Wrote: {OUT_CSV} rows={len(inv)}")
    print("Unique PDB IDs:", inv["pdb_id"].nunique() if len(inv) else 0)
    print("\nTop 10 UniProt accessions by #PDBs:")
    if len(inv):
        print(inv["uniprot_accession"].value_counts().head(10).to_string())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
