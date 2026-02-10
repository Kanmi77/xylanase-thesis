#!/usr/bin/env python3
"""
Fetch protein FASTA sequences from NCBI for RefSeq accessions in refseq_inventory.csv.

Input:
  data/curated/refseq_inventory.csv

Outputs:
  data/raw/refseq/refseq_proteins.fasta
  data/curated/refseq_fetch_status.csv

NCBI endpoint:
  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi

Notes:
- We batch IDs to reduce requests.
- If your VM/network blocks NCBI, status CSV will show failures.
"""

import os, time
import pandas as pd
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError

BASE = os.path.expanduser("~/xylanase-thesis")
INV = os.path.join(BASE, "data/curated/refseq_inventory.csv")
OUT_FASTA = os.path.join(BASE, "data/raw/refseq/refseq_proteins.fasta")
OUT_STATUS = os.path.join(BASE, "data/curated/refseq_fetch_status.csv")
LOG = os.path.join(BASE, "logs/07_fetch_refseq_fasta.log")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def fetch_fasta(batch_ids, timeout=180) -> tuple[bool, str]:
    params = {
        "db": "protein",
        "rettype": "fasta",
        "retmode": "text",
        "id": ",".join(batch_ids),
    }
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "xylanase-thesis/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            text = r.read().decode("utf-8", errors="replace")
        return True, text
    except (HTTPError, URLError) as e:
        return False, str(e)

def main() -> int:
    os.makedirs(os.path.dirname(OUT_FASTA), exist_ok=True)

    inv = pd.read_csv(INV, dtype=str).fillna("")
    ids = sorted(set(inv["refseq_acc"].tolist()))
    log(f"Unique RefSeq IDs: {len(ids)}")

    status_rows = []
    batch_size = 100  # safe
    ok_batches = 0
    fail_batches = 0

    # write FASTA incrementally
    with open(OUT_FASTA, "w", encoding="utf-8") as out_f:
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i+batch_size]
            ok, payload = fetch_fasta(batch)
            if ok and payload.strip().startswith(">"):
                out_f.write(payload)
                ok_batches += 1
                for rid in batch:
                    status_rows.append({"refseq_acc": rid, "status": "ok"})
            else:
                fail_batches += 1
                log(f"Batch failed at {i}..{i+len(batch)} : {payload}")
                for rid in batch:
                    status_rows.append({"refseq_acc": rid, "status": "fail"})
            time.sleep(0.35)  # be polite to NCBI

    pd.DataFrame(status_rows).to_csv(OUT_STATUS, index=False)
    log(f"Done. ok_batches={ok_batches} fail_batches={fail_batches}")
    log(f"Wrote: {OUT_FASTA}")
    log(f"Wrote: {OUT_STATUS}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
