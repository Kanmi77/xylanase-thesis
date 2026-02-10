#!/usr/bin/env python3
"""
Build a unified curated dataset:
- UniProt master (with organism_type + xref_cazy)
- RefSeq FASTA (downloaded) linked via refseq_inventory

Outputs:
  data/curated/xylanase_master_all_curated.csv

Rules:
- Keep only organism_type in {bacterial, fungal}
- Keep CAZy-derived GH label where available (GH10/GH11); else "unknown"
- Deduplicate by exact sequence (keep UniProt rows first when duplicates exist)
"""

import os, re
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")

UNIPROT_MASTER = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")
REFSEQ_INV = os.path.join(BASE, "data/curated/refseq_inventory.csv")
REFSEQ_FASTA = os.path.join(BASE, "data/raw/refseq/refseq_proteins.fasta")

OUT_ALL = os.path.join(BASE, "data/curated/xylanase_master_all_curated.csv")

GH10_RE = re.compile(r"\bGH10\b", re.IGNORECASE)
GH11_RE = re.compile(r"\bGH11\b", re.IGNORECASE)

def parse_fasta(path: str) -> dict:
    """Return dict: accession -> sequence for a FASTA with headers starting with >ACC ..."""
    seqs = {}
    cur_id = None
    cur = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id and cur:
                    seqs[cur_id] = "".join(cur).replace(" ", "").replace("\t", "")
                cur_id = line[1:].split()[0]
                cur = []
            else:
                cur.append(line)
        if cur_id and cur:
            seqs[cur_id] = "".join(cur).replace(" ", "").replace("\t", "")
    return seqs

def gh_from_cazy(xref_cazy: str) -> str:
    s = str(xref_cazy or "")
    if GH10_RE.search(s):
        return "GH10"
    if GH11_RE.search(s):
        return "GH11"
    return "unknown"

def main() -> int:
    os.makedirs(os.path.dirname(OUT_ALL), exist_ok=True)

    # ---------- UniProt master ----------
    u = pd.read_csv(UNIPROT_MASTER, dtype=str).fillna("")
    u["source"] = "uniprot"
    u["gh_family"] = u["xref_cazy"].apply(gh_from_cazy) if "xref_cazy" in u.columns else "unknown"

    # Keep fungi + bacteria only
    u = u[u["organism_type"].isin(["bacterial", "fungal"])].copy()

    u_out = pd.DataFrame({
        "source": u["source"],
        "primary_id": u["uniprot_accession"],
        "uniprot_accession": u["uniprot_accession"],
        "refseq_acc": "",
        "organism_type": u["organism_type"],
        "organism": u["organism"],
        "gh_family": u["gh_family"],
        "sequence": u["sequence"],
    })

    # ---------- RefSeq sequences ----------
    inv = pd.read_csv(REFSEQ_INV, dtype=str).fillna("")
    inv = inv[inv["organism_type"].isin(["bacterial", "fungal"])].copy()

    refseq_seqs = parse_fasta(REFSEQ_FASTA)

    # attach sequences (only those successfully fetched)
    inv["sequence"] = inv["refseq_acc"].map(refseq_seqs).fillna("")
    inv = inv[inv["sequence"].str.strip() != ""].copy()

    r_out = pd.DataFrame({
        "source": "refseq",
        "primary_id": inv["refseq_acc"],
        "uniprot_accession": inv["uniprot_accession"],
        "refseq_acc": inv["refseq_acc"],
        "organism_type": inv["organism_type"],
        "organism": inv["organism"],
        "gh_family": "unknown",  # GH label not provided by RefSeq download
        "sequence": inv["sequence"],
    })

    # ---------- Combine ----------
    all_df = pd.concat([u_out, r_out], ignore_index=True)

    # Deduplicate by sequence, keeping UniProt first (so CAZy GH labels are preserved where possible)
    all_df["source_rank"] = all_df["source"].map({"uniprot": 0, "refseq": 1}).fillna(9).astype(int)
    all_df = all_df.sort_values(["source_rank", "primary_id"]).drop_duplicates(subset=["sequence"], keep="first")
    all_df = all_df.drop(columns=["source_rank"])

    all_df.to_csv(OUT_ALL, index=False)

    print(f"Wrote: {OUT_ALL} rows={len(all_df)}")
    print("Counts by organism_type:")
    print(all_df["organism_type"].value_counts().to_string())
    print("\nCounts by source:")
    print(all_df["source"].value_counts().to_string())
    print("\nCounts by gh_family:")
    print(all_df["gh_family"].value_counts().to_string())

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
