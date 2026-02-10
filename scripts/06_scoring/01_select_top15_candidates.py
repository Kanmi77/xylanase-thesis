#!/usr/bin/env python3
import os
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
INP = os.path.join(BASE, "data/curated/xylanase_structured_subset_with_foldx_norm.csv")
OUT = os.path.join(BASE, "results/ranking/top15_candidates.csv")

# balanced quotas (adjust if you want)
QUOTAS = {
    ("bacterial", "GH10"): 4,
    ("bacterial", "GH11"): 4,
    ("fungal", "GH10"): 4,
    ("fungal", "GH11"): 3,
}

def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    df = pd.read_csv(INP, dtype=str).fillna("")
    # numeric
    df["foldx_energy_per_residue"] = pd.to_numeric(df["foldx_energy_per_residue"], errors="coerce")
    df["salt_bridge_count"] = pd.to_numeric(df["salt_bridge_count"], errors="coerce")
    df["disulfide_count"] = pd.to_numeric(df["disulfide_count"], errors="coerce")
    df["sasa_total"] = pd.to_numeric(df["sasa_total"], errors="coerce")
    df["chain_length"] = pd.to_numeric(df["chain_length"], errors="coerce")

    # keep only labeled GH and valid foldx
    df = df[df["gh_family"].isin(["GH10", "GH11"])]
    df = df[df["organism_type"].isin(["bacterial", "fungal"])]
    df = df[df["foldx_energy_per_residue"].notna()].copy()

    # optional: compactness proxy
    df["sasa_per_residue"] = df["sasa_total"] / df["chain_length"]

    # ranking: more negative is better; then more salt bridges; then disulfides; then lower sasa/res
    df = df.sort_values(
        by=["foldx_energy_per_residue", "salt_bridge_count", "disulfide_count", "sasa_per_residue"],
        ascending=[True, False, False, True]
    )

    picked = []
    for (org, gh), k in QUOTAS.items():
        sub = df[(df["organism_type"] == org) & (df["gh_family"] == gh)].head(k)
        picked.append(sub)

    top = pd.concat(picked, ignore_index=True)

    # final sort purely by stability after balancing
    top = top.sort_values(by=["foldx_energy_per_residue"], ascending=True)

    cols = [
        "pdb_id","chosen_chain","file_path","chain_length","organism_type","organism",
        "gh_family","xref_cazy","foldx_wt_total_energy","foldx_energy_per_residue",
        "salt_bridge_count","disulfide_count","sasa_total","sasa_per_residue","uniprot_accession","pdb_tag"
    ]
    top[cols].to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(top)}")
    print(top[["organism_type","gh_family"]].value_counts().to_string())

if __name__ == "__main__":
    main()
