#!/usr/bin/env python3
"""
Prepare FoldX-ready PDBs by extracting the chosen chain from each structure.

Input:
  data/curated/xylanase_structured_subset.csv
  data/raw/pdb/*.cif|*.pdb

Output:
  foldx/pdb_chains/<pdb_id>_<chain>.pdb
  data/curated/foldx_inputs.csv
"""

import os
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "data/curated/xylanase_structured_subset.csv")
OUT_DIR = os.path.join(BASE, "foldx/pdb_chains")
OUT_MAP = os.path.join(BASE, "data/curated/foldx_inputs.csv")
LOG = os.path.join(BASE, "logs/18A_prepare_foldx_chain_pdbs.log")

class ChainSelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return 1 if chain.id == self.chain_id else 0

def log(msg: str) -> None:
    print(msg)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

def load_structure(pdb_id: str, path: str):
    if path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, path)

def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(IN_CSV, dtype=str).fillna("")

    rows = []
    ok = 0
    fail = 0

    for _, r in df.iterrows():
        pdb_id = r["pdb_id"].upper()
        chain = r["chosen_chain"]
        path = r["file_path"]

        out_pdb = os.path.join(OUT_DIR, f"{pdb_id}_{chain}.pdb")

        if os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 0:
            rows.append({**r.to_dict(), "foldx_pdb": out_pdb})
            continue

        try:
            s = load_structure(pdb_id, path)
            io = PDBIO()
            io.set_structure(s)
            io.save(out_pdb, select=ChainSelect(chain))
            if os.path.getsize(out_pdb) > 0:
                ok += 1
                rows.append({**r.to_dict(), "foldx_pdb": out_pdb})
            else:
                fail += 1
                log(f"FAIL empty output: {pdb_id} chain {chain}")
        except Exception as e:
            fail += 1
            log(f"FAIL {pdb_id} chain {chain} -> {e}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT_MAP, index=False)

    log(f"Wrote: {OUT_MAP} rows={len(out)}")
    log(f"Prepared PDBs: ok={ok} fail={fail} (existing reused not counted)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
