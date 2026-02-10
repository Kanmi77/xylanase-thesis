#!/usr/bin/env python3
"""
Merge FoldX WT total energy into the structured subset.

Inputs:
  data/curated/xylanase_structured_subset.csv
  results/foldx/foldx4_wt_stability.csv

Output:
  data/curated/xylanase_structured_subset_with_foldx.csv
"""

import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
STRUCT = os.path.join(BASE, "data/curated/xylanase_structured_subset.csv")
FOLDX = os.path.join(BASE, "results/foldx/foldx4_wt_stability.csv")
OUT = os.path.join(BASE, "data/curated/xylanase_structured_subset_with_foldx.csv")

def main() -> int:
    s = pd.read_csv(STRUCT, dtype=str).fillna("")
    f = pd.read_csv(FOLDX, dtype=str).fillna("")

    # tag format matches your foldx file: <pdb_id>_<chain>_Repair
    s["pdb_tag"] = s["pdb_id"].str.upper() + "_" + s["chosen_chain"] + "_Repair"

    f["foldx_wt_total_energy"] = pd.to_numeric(f["foldx_wt_total_energy"], errors="coerce")

    m = s.merge(f[["pdb_tag", "foldx_wt_total_energy"]], on="pdb_tag", how="left")
    m.to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(m)} cols={m.shape[1]}")
    print("With FoldX energy:", int(m["foldx_wt_total_energy"].notna().sum()))
    print("Missing FoldX energy:", int(m["foldx_wt_total_energy"].isna().sum()))

    # list missing tags for traceability
    missing = m[m["foldx_wt_total_energy"].isna()][["pdb_id","chosen_chain","pdb_tag","file_path"]]
    if len(missing) > 0:
        miss_out = os.path.join(BASE, "results/foldx/foldx_missing_energy.csv")
        missing.to_csv(miss_out, index=False)
        print("Wrote missing list:", miss_out, "rows=", len(missing))

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
