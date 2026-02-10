#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "results/ranking/top15_candidates_with_docking.csv")
OUT_CSV = os.path.join(BASE, "results/ranking/top15_final_ranked.csv")

def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    mu = s.mean(skipna=True)
    sd = s.std(skipna=True)
    if sd == 0 or np.isnan(sd):
        return (s * 0.0)
    return (s - mu) / sd

def main():
    if not os.path.exists(IN_CSV):
        raise SystemExit(f"ERROR: Missing input: {IN_CSV}")

    df = pd.read_csv(IN_CSV, dtype=str).fillna("")

    # Numerics
    df["foldx_energy_per_residue"] = pd.to_numeric(df.get("foldx_energy_per_residue"), errors="coerce")
    df["salt_bridge_count"] = pd.to_numeric(df.get("salt_bridge_count"), errors="coerce")
    df["disulfide_count"] = pd.to_numeric(df.get("disulfide_count"), errors="coerce")
    df["vina_best_min"] = pd.to_numeric(df.get("vina_best_min"), errors="coerce")

    # Components:
    # FoldX: more negative is better -> multiply by -1 so bigger is better
    df["stability_component_raw"] = -1.0 * df["foldx_energy_per_residue"]

    # Docking: more negative is better -> multiply by -1 so bigger is better
    df["docking_component_raw"] = -1.0 * df["vina_best_min"]

    # Z-score normalize each component across the 15 candidates
    df["z_stability"] = zscore(df["stability_component_raw"])
    df["z_docking"] = zscore(df["docking_component_raw"])
    df["z_salt"] = zscore(df["salt_bridge_count"])
    df["z_disulfide"] = zscore(df["disulfide_count"])

    # Weights (simple + explainable)
    w_stab = 0.55
    w_dock = 0.20
    w_salt = 0.15
    w_dis = 0.10

    df["final_score"] = (
        w_stab * df["z_stability"] +
        w_dock * df["z_docking"] +
        w_salt * df["z_salt"] +
        w_dis  * df["z_disulfide"]
    )

    # Sort best first
    out = df.sort_values("final_score", ascending=False)

    # Keep a clean output view (but retain key columns)
    keep_cols = [
        "final_score",
        "pdb_id","chosen_chain","pdb_tag",
        "organism_type","gh_family","organism","uniprot_accession",
        "foldx_energy_per_residue","foldx_wt_total_energy",
        "salt_bridge_count","disulfide_count",
        "vina_best_xylobiose","vina_best_xylotriose","vina_best_min","vina_best_mean",
        "sasa_total","sasa_per_residue"
    ]
    keep_cols = [c for c in keep_cols if c in out.columns]
    out = out[keep_cols]

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    print(f"Wrote: {OUT_CSV} rows={len(out)}")
    print("\nTop 10 (preview):")
    print(out.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
