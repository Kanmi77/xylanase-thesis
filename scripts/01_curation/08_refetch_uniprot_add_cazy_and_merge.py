#!/usr/bin/env python3
"""
Refetch UniProt stream with CAZy column and merge into existing master CSV.

Inputs:
  data/raw/uniprot_xylanase.tsv
  data/curated/xylanase_master_uniprot.csv

Outputs:
  data/raw/uniprot_xylanase_with_cazy.tsv
  data/curated/xylanase_master_uniprot.csv  (updated with xref_cazy column)

Why:
- Your current TSV headers include only: Entry, Entry Name, Protein names, Gene Names, Organism,
  Organism (ID), Taxonomic lineage, Length, Sequence, PDB, RefSeq
- Proposal requires CAZy as a source, so we add CAZy xref.
"""

import os, time
import pandas as pd
import urllib.parse, urllib.request
from urllib.error import HTTPError

BASE = os.path.expanduser("~/xylanase-thesis")
OUT_TSV = os.path.join(BASE, "data/raw/uniprot_xylanase_with_cazy.tsv")
MASTER = os.path.join(BASE, "data/curated/xylanase_master_uniprot.csv")
LOG = os.path.join(BASE, "logs/08_refetch_uniprot_add_cazy_and_merge.log")

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def http_get(url: str, timeout: int = 180) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "xylanase-thesis/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        log(f"HTTPError {e.code}: {e.reason}")
        log(body[:1000])
        raise

def main() -> int:
    # minimal stable fields + CAZy
    query = "ec:3.2.1.8"
    fields = [
        "accession",
        "xref_cazy",
    ]
    params = {"query": query, "format": "tsv", "fields": ",".join(fields)}
    url = UNIPROT_STREAM + "?" + urllib.parse.urlencode(params)
    log(f"Fetching: {url}")

    data = http_get(url)
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "wb") as f:
        f.write(data)
    log(f"Wrote: {OUT_TSV} bytes={len(data)}")

    # Load and merge into master
    add = pd.read_csv(OUT_TSV, sep="\t", dtype=str).fillna("")
    # UniProt headers might come back as: Entry, CAZy (or similar)
    # detect accession column
    acc_col = "Entry" if "Entry" in add.columns else add.columns[0]
    cazy_col = [c for c in add.columns if c != acc_col][0]

    add = add.rename(columns={acc_col: "uniprot_accession", cazy_col: "xref_cazy"})
    master = pd.read_csv(MASTER, dtype=str).fillna("")
    merged = master.merge(add[["uniprot_accession", "xref_cazy"]], on="uniprot_accession", how="left")
    merged["xref_cazy"] = merged["xref_cazy"].fillna("").astype(str)

    merged.to_csv(MASTER, index=False)
    log("Merged CAZy into master CSV (updated in place).")
    log(f"Non-empty xref_cazy: {(merged['xref_cazy'].str.strip()!='').sum()}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
