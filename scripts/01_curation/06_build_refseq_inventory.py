#!/usr/bin/env python3
"""
Parse RefSeq accessions from the curated UniProt master CSV.

Input:
  data/curated/xylanase_master_uniprot.csv  (column: xref_refseq)

Output:
  data/curated/refseq_inventory.csv
    uniprot_accession, organism_type, organism, refseq_acc
"""

import os, re
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")
OUT_CSV = os.path.join(BASE, "data/curated/refseq_inventory.csv")

# RefSeq protein accession patterns (common)
REFSEQ_RE = re.compile(r"\b(?:NP|XP|WP|YP|AP)_\d+(?:\.\d+)?\b")

def main() -> int:
    df = pd.read_csv(IN_CSV, dtype=str).fillna("")
    rows = []

    for _, r in df.iterrows():
        acc = r["uniprot_accession"].strip()
        ref = str(r.get("xref_refseq", ""))
        ids = sorted(set(REFSEQ_RE.findall(ref)))
        for rid in ids:
            rows.append({
                "uniprot_accession": acc,
                "organism_type": r.get("organism_type", ""),
                "organism": r.get("organism", ""),
                "refseq_acc": rid
            })

    inv = pd.DataFrame(rows).drop_duplicates()
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    inv.to_csv(OUT_CSV, index=False)

    print(f"Wrote: {OUT_CSV} rows={len(inv)} unique_refseq={inv['refseq_acc'].nunique() if len(inv) else 0}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
