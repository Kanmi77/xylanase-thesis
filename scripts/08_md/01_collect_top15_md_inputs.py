#!/usr/bin/env python3
import os, shutil
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")

TOP15 = os.path.join(BASE, "results/ranking/top15_final_ranked.csv")
DOCK_MANIFEST = os.path.join(BASE, "docking/manifests/dock_manifest_top15.csv")
IN_DIR = os.path.join(BASE, "docking/proteins")
OUT_DIR = os.path.join(BASE, "md/inputs")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(DOCK_MANIFEST):
        raise SystemExit(f"Missing: {DOCK_MANIFEST}")
    if not os.path.exists(TOP15):
        raise SystemExit(f"Missing: {TOP15}")

    df = pd.read_csv(DOCK_MANIFEST, dtype=str).fillna("")
    top = pd.read_csv(TOP15, dtype=str).fillna("")
    top_tags = set(top["pdb_tag"].tolist())

    df = df[df["pdb_tag"].isin(top_tags)].copy()
    if df.empty:
        raise SystemExit("No matching pdb_tags between top15_final_ranked.csv and dock_manifest_top15.csv")

    copied = 0
    for _, r in df.iterrows():
        tag = r["pdb_tag"]
        src = os.path.join(IN_DIR, f"{tag}.protein_clean.fixed.pdb")
        if not os.path.exists(src):
            raise SystemExit(f"Missing input PDB for {tag}: {src}")
        dst = os.path.join(OUT_DIR, f"{tag}.pdb")
        shutil.copyfile(src, dst)
        copied += 1

    print(f"Copied {copied} PDBs to: {OUT_DIR}")
    print("Example:", sorted(os.listdir(OUT_DIR))[:5])

if __name__ == "__main__":
    main()
