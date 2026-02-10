#!/usr/bin/env python3
import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")

TOP15_IN = os.path.join(BASE, "results/ranking/top15_candidates.csv")
DOCK_SUM = os.path.join(BASE, "results/docking/vina_top15_summary.csv")
OUT = os.path.join(BASE, "results/ranking/top15_candidates_with_docking.csv")

def main():
    if not os.path.exists(TOP15_IN):
        raise SystemExit(f"ERROR: Missing top15 file: {TOP15_IN}")
    if not os.path.exists(DOCK_SUM):
        raise SystemExit(f"ERROR: Missing docking summary: {DOCK_SUM}")

    top = pd.read_csv(TOP15_IN, dtype=str).fillna("")
    dock = pd.read_csv(DOCK_SUM, dtype=str).fillna("")

    # Convert numeric fields that matter
    for c in ["foldx_energy_per_residue", "salt_bridge_count", "disulfide_count", "sasa_per_residue"]:
        if c in top.columns:
            top[c] = pd.to_numeric(top[c], errors="coerce")

    for c in ["vina_best_xylobiose", "vina_best_xylotriose", "vina_best_min", "vina_best_mean"]:
        if c in dock.columns:
            dock[c] = pd.to_numeric(dock[c], errors="coerce")

    # Merge on pdb_tag (unique identifier for receptor structure)
    merged = top.merge(
        dock[["pdb_tag","vina_best_xylobiose","vina_best_xylotriose","vina_best_min","vina_best_mean"]],
        on="pdb_tag",
        how="left"
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    merged.to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(merged)}")
    print("Docking rows missing:", merged["vina_best_min"].isna().sum())

if __name__ == "__main__":
    main()
