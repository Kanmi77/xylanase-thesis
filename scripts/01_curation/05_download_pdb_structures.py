#!/usr/bin/env python3
"""
Download PDB structures for IDs in pdb_inventory.csv.

Downloads:
  - mmCIF (preferred): https://files.rcsb.org/download/<PDB>.cif
  - fallback PDB:      https://files.rcsb.org/download/<PDB>.pdb

Output directory:
  data/raw/pdb/

Log:
  logs/05_download_pdb_structures.log
"""

import os, time
import pandas as pd
import urllib.request
from urllib.error import HTTPError, URLError

BASE = os.path.expanduser("~/xylanase-thesis")
INV = os.path.join(BASE, "data/curated/pdb_inventory.csv")
OUTDIR = os.path.join(BASE, "data/raw/pdb")
LOG = os.path.join(BASE, "logs/05_download_pdb_structures.log")

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def fetch(url: str, outpath: str, timeout: int = 180) -> bool:
    req = urllib.request.Request(url, headers={"User-Agent": "xylanase-thesis/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
        with open(outpath, "wb") as f:
            f.write(data)
        return True
    except (HTTPError, URLError) as e:
        log(f"FAILED {url} -> {e}")
        return False

def main() -> int:
    os.makedirs(OUTDIR, exist_ok=True)

    inv = pd.read_csv(INV, dtype=str).fillna("")
    pdb_ids = sorted(set(inv["pdb_id"].str.upper().tolist()))

    log(f"Total unique PDB IDs: {len(pdb_ids)}")
    ok, fail, skip = 0, 0, 0

    for pid in pdb_ids:
        cif_path = os.path.join(OUTDIR, f"{pid}.cif")
        pdb_path = os.path.join(OUTDIR, f"{pid}.pdb")

        if os.path.exists(cif_path) or os.path.exists(pdb_path):
            skip += 1
            continue

        cif_url = f"https://files.rcsb.org/download/{pid}.cif"
        pdb_url = f"https://files.rcsb.org/download/{pid}.pdb"

        if fetch(cif_url, cif_path):
            ok += 1
        elif fetch(pdb_url, pdb_path):
            ok += 1
        else:
            fail += 1

        time.sleep(0.15)

    log(f"Done. ok={ok} fail={fail} skip_existing={skip}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
