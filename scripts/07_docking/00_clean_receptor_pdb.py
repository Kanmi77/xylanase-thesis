#!/usr/bin/env python3
"""
Robust receptor cleaner for Meeko/RDKit.

Problem:
- mk_prepare_receptor.py fails if certain metals (e.g., Cd) are present.
- Some structures may label metals as ATOM (not HETATM), so "keep ATOM only" is not sufficient.

Solution:
- Keep protein ATOM lines BUT remove any atom whose element matches a banned set.
- Remove all HETATM by default (waters, ligands, ions).
- If element column is blank, infer element from atom name.
- Also remove waters (HOH/WAT) even if mislabelled.

Output is a PDB safe for mk_prepare_receptor.py.

Usage:
  00_clean_receptor_pdb.py input.pdb output.protein_clean.pdb
"""

import sys
import re

# Metals/ions that commonly break bond perception or are unwanted for basic docking receptor prep.
BANNED_ELEMENTS = {
    "CD", "HG", "ZN", "FE", "MN", "CU", "CO", "NI",
    "CA", "MG", "NA", "K", "CL", "BR", "I",
    "SR", "CS", "BA", "PB", "SN", "AG", "AU"
}
# NOTE: "CA" element means calcium ion. Protein alpha-carbon atom name "CA" is carbon (C), not element CA.
# Our inference avoids labeling protein CA as calcium.

AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"
}

TWO_LETTER_ELEMENTS = {
    "CL","BR","NA","MG","ZN","FE","MN","CO","NI","CU","CD","HG","SR","CS","BA",
    "AL","LI","RB","AG","AU","PB","SN","SE","AS","CR","TI","CA","SI"
}

def infer_element(atom_name: str, res_name: str) -> str:
    a = atom_name.strip()
    r = res_name.strip().upper()

    if not a:
        return ""

    # remove leading digits like 1HG1
    a2 = re.sub(r"^[0-9]+", "", a).upper()
    if not a2:
        return ""

    # protein alpha-carbon
    if a2 == "CA" and (r in AA3 or r == ""):
        return "C"

    if len(a2) >= 2 and a2[:2] in TWO_LETTER_ELEMENTS:
        return a2[:2].upper()

    return a2[0].upper()

def get_element_from_line(line: str) -> str:
    # element col 77-78 (0-based 76:78)
    elem = (line[76:78] if len(line) >= 78 else "").strip().upper()
    if elem:
        return elem
    atom_name = line[12:16] if len(line) >= 16 else ""
    res_name  = line[17:20] if len(line) >= 20 else ""
    return infer_element(atom_name, res_name)

def is_water(line: str) -> bool:
    res = (line[17:20] if len(line) >= 20 else "").strip().upper()
    return res in {"HOH","WAT","H2O","DOD"}

def clean(inp: str, outp: str) -> None:
    kept = 0
    dropped_hetatm = 0
    dropped_banned = 0
    dropped_water = 0
    total_atoms = 0

    with open(inp, "r", encoding="utf-8", errors="ignore") as fin, \
         open(outp, "w", encoding="utf-8") as fout:

        for line in fin:
            if line.startswith(("ATOM  ", "HETATM")):
                total_atoms += 1

                # Drop all HETATM outright (standard docking receptor prep)
                if line.startswith("HETATM"):
                    dropped_hetatm += 1
                    continue

                # Drop water even if mislabelled
                if is_water(line):
                    dropped_water += 1
                    continue

                # Drop banned elements (e.g., Cd)
                elem = get_element_from_line(line)
                if elem in BANNED_ELEMENTS:
                    dropped_banned += 1
                    continue

                # Otherwise keep ATOM
                fout.write(line)
                kept += 1

            elif line.startswith(("TER", "END")):
                fout.write(line)
            else:
                # ignore headers/remarks
                pass

    print(f"Cleaned receptor: {inp} -> {outp}")
    print(f"  total ATOM/HETATM lines seen: {total_atoms}")
    print(f"  kept ATOM lines: {kept}")
    print(f"  dropped HETATM lines: {dropped_hetatm}")
    print(f"  dropped water lines: {dropped_water}")
    print(f"  dropped banned-element lines: {dropped_banned}")

def main():
    if len(sys.argv) != 3:
        print("Usage: 00_clean_receptor_pdb.py input.pdb output.protein_clean.pdb", file=sys.stderr)
        return 2
    clean(sys.argv[1], sys.argv[2])
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
