#!/usr/bin/env python3
import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")

IN_CSV = os.path.join(BASE, "results/docking/vina_top15_scores.csv")
OUT_CSV = os.path.join(BASE, "results/docking/vina_top15_summary.csv")

def main():
    if not os.path.exists(IN_CSV):
        raise SystemExit(f"ERROR: Missing input: {IN_CSV}")

    df = pd.read_csv(IN_CSV, dtype=str).fillna("")
    if df.empty:
        raise SystemExit("ERROR: vina_top15_scores.csv is empty")

    # Convert affinity to float
    df["best_affinity_kcal_mol"] = df["best_affinity_kcal_mol"].astype(float)

    # Keep only mode-1 best scores (your parser already did this)
    # Pivot by ligand so we get 1 row per pdb_tag
    piv = df.pivot_table(
        index=["pdb_tag","pdb_id","chosen_chain","organism_type","organism","uniprot_accession","gh_family"],
        columns="ligand",
        values="best_affinity_kcal_mol",
        aggfunc="min"  # in case duplicates
    ).reset_index()

    # Ensure expected ligand columns exist
    for lig in ["xylobiose", "xylotriose"]:
        if lig not in piv.columns:
            piv[lig] = pd.NA

    piv = piv.rename(columns={
        "xylobiose": "vina_best_xylobiose",
        "xylotriose": "vina_best_xylotriose"
    })

    # Combined metrics
    piv["vina_best_min"] = piv[["vina_best_xylobiose", "vina_best_xylotriose"]].min(axis=1, skipna=True)
    piv["vina_best_mean"] = piv[["vina_best_xylobiose", "vina_best_xylotriose"]].mean(axis=1, skipna=True)

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    piv.to_csv(OUT_CSV, index=False)

    print(f"Wrote: {OUT_CSV} rows={len(piv)}")
    print("\nSanity checks:")
    print("Unique pdb_tag:", piv["pdb_tag"].nunique())
    print("Missing xylobiose:", piv["vina_best_xylobiose"].isna().sum())
    print("Missing xylotriose:", piv["vina_best_xylotriose"].isna().sum())

if __name__ == "__main__":
    main()
