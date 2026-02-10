#!/usr/bin/env python3
"""
Export 4 FASTA files from xylanase_master_all_curated.csv:

- fungal_GH10.fasta
- fungal_GH11.fasta
- bacterial_GH10.fasta
- bacterial_GH11.fasta

Only includes rows with gh_family in {GH10, GH11}.
"""

import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "data/curated/xylanase_master_all_curated.csv")
OUT_DIR = os.path.join(BASE, "data/curated/splits_fasta")

def write_fasta(df: pd.DataFrame, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for _, r in df.iterrows():
            pid = str(r["primary_id"]).strip()
            src = str(r["source"]).strip()
            org = str(r["organism"]).strip()
            seq = str(r["sequence"]).strip()
            if not pid or not seq:
                continue
            f.write(f">{pid} | {src} | {org}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(IN_CSV, dtype=str).fillna("")

    df = df[df["gh_family"].isin(["GH10", "GH11"])].copy()

    for org_type in ["fungal", "bacterial"]:
        for gh in ["GH10", "GH11"]:
            sub = df[(df["organism_type"] == org_type) & (df["gh_family"] == gh)].copy()
            outpath = os.path.join(OUT_DIR, f"{org_type}_{gh}.fasta")
            write_fasta(sub, outpath)
            print(f"Wrote: {outpath} entries={len(sub)}")

    # Summary
    print("\nSummary counts:")
    print(df.groupby(["organism_type","gh_family"]).size().to_string())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
