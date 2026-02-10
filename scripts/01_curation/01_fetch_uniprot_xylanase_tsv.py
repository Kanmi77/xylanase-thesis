#!/usr/bin/env python3
"""
Fetch UniProtKB entries for xylanase (EC 3.2.1.8) in TSV format using UniProt REST streaming API.

This version uses a minimal set of UniProt fields that are reliably supported by the TSV streamer.
We can expand fields after this works.

Outputs:
  data/raw/uniprot_xylanase.tsv
Logs:
  logs/01_fetch_uniprot_xylanase_tsv.log
"""
import os, time
import urllib.parse
import urllib.request
from urllib.error import HTTPError

BASE = os.path.expanduser("~/xylanase-thesis")
OUT_TSV = os.path.join(BASE, "data/raw/uniprot_xylanase.tsv")
LOG = os.path.join(BASE, "logs/01_fetch_uniprot_xylanase_tsv.log")

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
        # print UniProt's error body (this is what we need for debugging)
        body = e.read().decode("utf-8", errors="replace")
        log(f"HTTPError {e.code}: {e.reason}")
        log("---- UniProt error body (first 1000 chars) ----")
        log(body[:1000])
        raise

def main() -> int:
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

    # Keep query simple and URL-safe
    query = "ec:3.2.1.8"

    # SAFE field set (expand later)
    fields = [
        "accession",
        "id",
        "protein_name",
        "gene_names",
        "organism_name",
        "organism_id",
        "lineage",
        "length",
        "sequence",
        "xref_pdb",
        "xref_refseq",
    ]

    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
    }

    url = UNIPROT_STREAM + "?" + urllib.parse.urlencode(params)
    log(f"Fetching UniProt stream: query=({query})")
    log(f"URL: {url}")

    data = http_get(url)
    with open(OUT_TSV, "wb") as f:
        f.write(data)

    log(f"Wrote: {OUT_TSV} ({len(data)} bytes)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
