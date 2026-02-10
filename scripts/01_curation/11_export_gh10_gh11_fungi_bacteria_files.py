#!/usr/bin/env python3
"""
Export 4 splits (fungal/bacterial x GH10/GH11) from UniProt master dataset.

Inputs:
  data/curated/xylanase_master_uniprot.csv

Outputs (8 files):
  data/curated/splits_gh/
    fungal_GH10.csv
    fungal_GH11.csv
    bacterial_GH10.csv
    bacterial_GH11.csv
    fungal_GH10.fasta
    fungal_GH11.fasta
    bacterial_GH10.fasta
    bacterial_GH11.fasta

GH family is derived from xref_cazy field (must contain GH10 or GH11).
"""

import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")
OUT_DIR = os.path.join(BASE, "data/curated/splits_gh")

def write_fasta(df: pd.DataFrame, outpath: str) -> None:
    with open(outpath, "w", encoding="utf-8") as f:
        for _, r in df.iterrows():
            acc = str(r["uniprot_accession"]).strip()
            org = str(r["organism"]).strip()
            typ = str(r["organism_type"]).strip()
            gh = str(r["gh_family"]).strip()
            seq = str(r["sequence"]).strip()
            if not acc or not seq:
                continue
            f.write(f">{acc} | {typ} | {gh} | {org}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)

    df = pd.read_csv(IN_CSV, dtype=str).fillna("")
    if "xref_cazy" not in df.columns:
        raise SystemExit("ERROR: xref_cazy column missing in UniProt master.")

    df["gh_family"] = "unknown"
    df.loc[df["xref_cazy"].str.contains("GH10", na=False), "gh_family"] = "GH10"
    df.loc[df["xref_cazy"].str.contains("GH11", na=False), "gh_family"] = "GH11"

    df = df[df["organism_type"].isin(["fungal","bacterial"])].copy()
    df = df[df["gh_family"].isin(["GH10","GH11"])].copy()

    # Keep only core columns for CSV outputs
    core_cols = [
        "uniprot_accession","uniprot_id","protein_name","gene_names",
        "organism_type","organism","taxonomy_id","length","gh_family",
        "xref_cazy","xref_pdb","xref_refseq","sequence"
    ]
    core_cols = [c for c in core_cols if c in df.columns]

    for org_type in ["fungal","bacterial"]:
        for gh in ["GH10","GH11"]:
            sub = df[(df["organism_type"]==org_type) & (df["gh_family"]==gh)].copy()
            csv_path = os.path.join(OUT_DIR, f"{org_type}_{gh}.csv")
            fa_path  = os.path.join(OUT_DIR, f"{org_type}_{gh}.fasta")

            sub[core_cols].to_csv(csv_path, index=False)
            write_fasta(sub, fa_path)

            print(f"Wrote: {csv_path} rows={len(sub)}")
            print(f"Wrote: {fa_path} entries={len(sub)}")

    print("\nSummary counts:")
    print(df.groupby(["organism_type","gh_family"]).size().to_string())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
