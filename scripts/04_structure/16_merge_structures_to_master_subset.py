#!/usr/bin/env python3
"""
Merge PDB manifest + structural features + UniProt mapping into a structured subset table.

Inputs:
  results/structures/structure_manifest.csv
  results/structures/structural_features.csv
  data/curated/pdb_inventory.csv
  data/curated/xylanase_master_uniprot.csv

Output:
  data/curated/xylanase_structured_subset.csv

Notes:
- One PDB can map to multiple UniProt accessions; we keep one row per (pdb_id, uniprot_accession).
"""

import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")

MANIFEST = os.path.join(BASE, "results/structures/structure_manifest.csv")
FEATURES = os.path.join(BASE, "results/structures/structural_features.csv")
PDB_INV  = os.path.join(BASE, "data/curated/pdb_inventory.csv")
MASTER   = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")

OUT = os.path.join(BASE, "data/curated/xylanase_structured_subset.csv")

def main() -> int:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    man = pd.read_csv(MANIFEST, dtype=str).fillna("")
    feat = pd.read_csv(FEATURES, dtype=str).fillna("")
    inv = pd.read_csv(PDB_INV, dtype=str).fillna("")
    master = pd.read_csv(MASTER, dtype=str).fillna("")

    # Normalize IDs
    man["pdb_id"] = man["pdb_id"].str.upper()
    feat["pdb_id"] = feat["pdb_id"].str.upper()
    inv["pdb_id"] = inv["pdb_id"].str.upper()

    # merge manifest + features
    sf = man.merge(feat, on=["pdb_id", "chosen_chain"], how="left")

    # map pdb -> uniprot
    sf = sf.merge(inv[["pdb_id", "uniprot_accession"]].drop_duplicates(),
                  on="pdb_id", how="left")

    # bring UniProt metadata
    keep_cols = ["uniprot_accession","organism_type","organism","xref_cazy"]
    master_small = master[keep_cols].drop_duplicates()
    sf = sf.merge(master_small, on="uniprot_accession", how="left")

    # GH family from CAZy xref
    sf["gh_family"] = "unknown"
    sf.loc[sf["xref_cazy"].str.contains("GH10", na=False), "gh_family"] = "GH10"
    sf.loc[sf["xref_cazy"].str.contains("GH11", na=False), "gh_family"] = "GH11"

    # tidy numeric columns (keep as numeric-friendly)
    num_cols = ["chain_length","hbond_proxy_count","salt_bridge_count","disulfide_count","sasa_total"]
    for c in num_cols:
        if c in sf.columns:
            sf[c] = pd.to_numeric(sf[c], errors="coerce")

    sf.to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(sf)} cols={len(sf.columns)}")
    print("\nCounts by organism_type:")
    print(sf["organism_type"].value_counts(dropna=False).to_string())
    print("\nCounts by gh_family:")
    print(sf["gh_family"].value_counts(dropna=False).to_string())
    print("\nUnique PDB IDs:", sf["pdb_id"].nunique())
    print("Unique UniProt accessions:", sf["uniprot_accession"].nunique())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
