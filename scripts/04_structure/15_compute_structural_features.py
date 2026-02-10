#!/usr/bin/env python3
"""
Compute structural features from structures:
- salt bridges (Asp/Glu vs Lys/Arg/His, cutoff 4.0 Å)
- disulfide bonds (Cys SG-SG cutoff 2.2 Å)
- H-bond proxy (N/O distance cutoff 3.5 Å)
- SASA total (Shrake-Rupley)

Inputs:
  results/structures/structure_manifest.csv

Outputs:
  results/structures/structural_features.csv
"""

import os, math
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

BASE = os.path.expanduser("~/xylanase-thesis")
MANIFEST = os.path.join(BASE, "results/structures/structure_manifest.csv")
OUT = os.path.join(BASE, "results/structures/structural_features.csv")

ACIDIC = {"ASP","GLU"}
BASIC = {"LYS","ARG","HIS"}

def dist(a, b) -> float:
    dx = a[0]-b[0]; dy=a[1]-b[1]; dz=a[2]-b[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def load_model(pid: str, path: str):
    parser = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    s = parser.get_structure(pid, path)
    return next(s.get_models())

def main() -> int:
    man = pd.read_csv(MANIFEST, dtype=str).fillna("")
    rows = []
    sr = ShrakeRupley()

    for _, r in man.iterrows():
        pid = r["pdb_id"]
        path = r["file_path"]
        chain_id = r["chosen_chain"]

        if not chain_id:
            continue

        model = load_model(pid, path)
        if chain_id not in model:
            continue

        chain = model[chain_id]

        # --- H-bond proxy: N/O atom pairs within 3.5 Å (exclude same residue) ---
        atoms = [a for a in chain.get_atoms() if a.element in ("N","O")]
        hb = 0
        for i in range(len(atoms)):
            ai = atoms[i]
            ri = ai.get_parent()
            for j in range(i+1, len(atoms)):
                aj = atoms[j]
                rj = aj.get_parent()
                if ri == rj:
                    continue
                if dist(ai.coord, aj.coord) <= 3.5:
                    hb += 1

        # --- Salt bridges: acidic O to basic N within 4.0 Å ---
        acidic_atoms = []
        basic_atoms = []
        for res in chain.get_residues():
            if res.id[0] != " ":
                continue
            if res.resname in ACIDIC:
                for an in res.get_atoms():
                    if an.element == "O":
                        acidic_atoms.append(an)
            elif res.resname in BASIC:
                for an in res.get_atoms():
                    if an.element == "N":
                        basic_atoms.append(an)

        sb = 0
        for a in acidic_atoms:
            for b in basic_atoms:
                if dist(a.coord, b.coord) <= 4.0:
                    sb += 1

        # --- Disulfides: CYS SG-SG within 2.2 Å ---
        sg = []
        for res in chain.get_residues():
            if res.id[0] != " ":
                continue
            if res.resname == "CYS" and "SG" in res:
                sg.append(res["SG"])
        ds = 0
        for i in range(len(sg)):
            for j in range(i+1, len(sg)):
                if dist(sg[i].coord, sg[j].coord) <= 2.2:
                    ds += 1

        # --- SASA total ---
        sasa_total = float("nan")
        try:
            sr.compute(chain, level="R")
            sasa_total = sum(getattr(res, "sasa", 0.0) for res in chain.get_residues() if res.id[0] == " ")
        except Exception:
            pass

        rows.append({
            "pdb_id": pid,
            "chosen_chain": chain_id,
            "hbond_proxy_count": hb,
            "salt_bridge_count": sb,
            "disulfide_count": ds,
            "sasa_total": sasa_total,
        })

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, index=False)

    print(f"Wrote: {OUT} rows={len(df)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
