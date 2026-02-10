#!/usr/bin/env python3
"""
Add FoldX energy per residue to the merged structured subset.

Input:
  data/curated/xylanase_structured_subset_with_foldx.csv

Output:
  data/curated/xylanase_structured_subset_with_foldx_norm.csv
"""

import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
INP = os.path.join(BASE, "data/curated/xylanase_structured_subset_with_foldx.csv")
OUT = os.path.join(BASE, "data/curated/xylanase_structured_subset_with_foldx_norm.csv")

def main() -> int:
    df = pd.read_csv(INP, dtype=str).fillna("")

    df["chain_length_num"] = pd.to_numeric(df["chain_length"], errors="coerce")
    df["foldx_wt_total_energy_num"] = pd.to_numeric(df["foldx_wt_total_energy"], errors="coerce")

    df["foldx_energy_per_residue"] = df["foldx_wt_total_energy_num"] / df["chain_length_num"]

    # drop helper columns
    df.drop(columns=["chain_length_num", "foldx_wt_total_energy_num"], inplace=True)

    df.to_csv(OUT, index=False)
    print(f"Wrote: {OUT} rows={len(df)}")
    print("Non-null energy_per_residue:", df["foldx_energy_per_residue"].notna().sum())
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
