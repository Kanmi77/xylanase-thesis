#!/usr/bin/env python3
"""
Fix missing/blank element symbols in PDB files by populating columns 77-78.

Meeko/RDKit requires valid element symbols. Many PDB writers leave element blank.
We infer element from atom name (columns 13-16).

Input:  PDB file path
Output: fixed PDB file path
"""

import sys
import re

# Common PDB element inference rules:
# - Atom name like "CA" could be Calcium OR Carbon-alpha. In proteins, " CA " is carbon-alpha and element is C.
# - If atom name starts with a digit (e.g., 1HG1), element is next letter(s).
# - Halogens: CL, BR
# - Metals: FE, MG, ZN, MN, CA (but protein " CA " is C)
#
# We use conservative protein-centric rules.

TWO_LETTER = {"CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CO", "NI", "CU", "CD", "HG", "SR", "CS", "BA", "AL", "LI", "RB", "AG", "AU", "PB", "SN", "SE", "AS", "CR", "TI", "CA"}

def infer_element(atom_name: str, res_name: str) -> str:
    """
    Infer element symbol from atom name and residue name.
    """
    a = atom_name.strip()

    if not a:
        return ""

    # PDB atom name may have leading digit: 1HG, 2HB etc.
    a2 = re.sub(r"^[0-9]+", "", a).upper()

    # Special case: protein alpha carbon is atom name "CA" in residues like ALA, GLY, etc.
    # That should be element C (not calcium).
    if a2 == "CA":
        # If residue looks like a standard amino acid or common protein residue, treat as carbon.
        # Residue name is columns 18-20.
        r = res_name.strip().upper()
        if r and len(r) == 3:
            return "C"

    # If first 2 letters form known element like CL, BR, FE, ZN...
    if len(a2) >= 2 and a2[:2] in TWO_LETTER:
        return a2[:2].title()

    # Otherwise element is first letter
    return a2[0].title()

def fix_file(inp: str, outp: str) -> int:
    changed = 0
    total = 0

    with open(inp, "r", encoding="utf-8", errors="ignore") as fin, open(outp, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM")):
                total += 1

                # Atom name is cols 13-16 (0-based 12:16)
                atom_name = line[12:16]
                res_name = line[17:20] if len(line) >= 20 else ""

                # Element is cols 77-78 (0-based 76:78)
                elem = line[76:78] if len(line) >= 78 else ""
                elem_stripped = elem.strip()

                if elem_stripped == "":
                    e = infer_element(atom_name, res_name)
                    # ensure line length at least 78
                    l = line.rstrip("\n")
                    if len(l) < 78:
                        l = l.ljust(78)
                    # write element right-aligned in 2 columns
                    l = l[:76] + f"{e:>2}" + l[78:]
                    line = l + "\n"
                    changed += 1

            fout.write(line)

    print(f"Fixed PDB elements: {inp} -> {outp}  (ATOM/HETATM={total}, changed={changed})")
    return 0

def main():
    if len(sys.argv) != 3:
        print("Usage: 00_fix_pdb_elements.py input.pdb output_fixed.pdb", file=sys.stderr)
        return 2
    return fix_file(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    raise SystemExit(main())
