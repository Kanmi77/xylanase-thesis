#!/usr/bin/env python3
"""
Curate UniProt TSV into a master CSV using the exact headers present in your file.

Input:
  data/raw/uniprot_xylanase.tsv  (headers observed)
    ['Entry','Entry Name','Protein names','Gene Names','Organism','Organism (ID)',
     'Taxonomic lineage','Length','Sequence','PDB','RefSeq']

Output:
  data/curated/xylanase_master_uniprot.csv
"""

import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_TSV = os.path.join(BASE, "data/raw/uniprot_xylanase.tsv")
OUT_CSV = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")

def is_missing(v: str) -> bool:
    if not isinstance(v, str):
        return True
    s = v.strip()
    return s == "" or s == "-" or s.lower() == "nan"

def classify_domain(lineage: str) -> str:
    if is_missing(lineage):
        return "unknown"
    if "Bacteria" in lineage:
        return "bacterial"
    if "Eukaryota" in lineage:
        if "Fungi" in lineage:
            return "fungal"
        return "eukaryote_non_fungal"
    if "Archaea" in lineage:
        return "archaea"
    return "unknown"

def has_pdb(v: str) -> int:
    return 0 if is_missing(v) else 1

def main() -> int:
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

    df = pd.read_csv(IN_TSV, sep="\t", dtype=str).fillna("")

    required = ["Entry", "Organism", "Taxonomic lineage", "Sequence"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: Missing required columns {missing}. Found: {df.columns.tolist()}")

    out = pd.DataFrame({
        "uniprot_accession": df["Entry"],
        "uniprot_id": df["Entry Name"] if "Entry Name" in df.columns else "",
        "protein_name": df["Protein names"] if "Protein names" in df.columns else "",
        "gene_names": df["Gene Names"] if "Gene Names" in df.columns else "",
        "organism": df["Organism"],
        "taxonomy_id": df["Organism (ID)"] if "Organism (ID)" in df.columns else "",
        "lineage": df["Taxonomic lineage"],
        "length": df["Length"] if "Length" in df.columns else "",
        "sequence": df["Sequence"],
        "xref_pdb": df["PDB"] if "PDB" in df.columns else "",
        "xref_refseq": df["RefSeq"] if "RefSeq" in df.columns else "",
    })

    # clean sequences
    out["sequence"] = out["sequence"].astype(str).str.replace(r"\s+", "", regex=True)

    # derived columns
    out["organism_type"] = out["lineage"].apply(classify_domain)
    out["has_pdb"] = out["xref_pdb"].apply(has_pdb)

    # drop empties
    out = out[(~out["uniprot_accession"].apply(is_missing)) & (~out["sequence"].apply(is_missing))]

    out.to_csv(OUT_CSV, index=False)

    print(f"Wrote: {OUT_CSV} rows={len(out)} cols={len(out.columns)}")
    print("\nOrganism_type counts:")
    print(out["organism_type"].value_counts().to_string())
    print("\nHas PDB counts:")
    print(out["has_pdb"].value_counts().to_string())
    print("\nExample PDB column values (first 15):")
    print(out["xref_pdb"].head(15).to_string(index=False))

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
