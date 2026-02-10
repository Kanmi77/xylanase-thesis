#!/usr/bin/env python3
"""
Sanitize a protein PDB so GROMACS pdb2gmx doesn't fail on incomplete residues.

Strategy (batch-safe for screening):
- Keep only ATOM records (protein)
- Drop altlocs other than ' ' or 'A'
- For each residue:
    - If residue is a standard amino acid but is missing required backbone atoms (N, CA, C, O): drop residue
    - If residue is standard AA but appears to have incomplete/odd sidechain:
        -> convert residue to GLY and keep only backbone atoms (N, CA, C, O, OXT)
This prevents errors like "ARG missing CD", etc.
"""

import argparse
from collections import defaultdict

BACKBONE = {"N", "CA", "C", "O", "OXT"}

# Standard amino acids in PDB 3-letter codes (common set)
STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}

def parse_atom_name(line: str) -> str:
    return line[12:16].strip()

def parse_altloc(line: str) -> str:
    return line[16:17]  # altLoc column

def parse_resname(line: str) -> str:
    return line[17:20].strip()

def parse_chain(line: str) -> str:
    return line[21:22]

def parse_resseq(line: str) -> int:
    # residue sequence number columns 23-26 (1-based indexing in docs)
    return int(line[22:26].strip())

def parse_icode(line: str) -> str:
    return line[26:27].strip()

def residue_key(line: str):
    return (parse_chain(line), parse_resseq(line), parse_icode(line))

def sanitize(in_pdb: str, out_pdb: str, keep_altloc=("","A")):
    with open(in_pdb, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # collect residues
    residues = defaultdict(list)
    header = []
    for ln in lines:
        rec = ln[0:6].strip()
        if rec in {"HEADER","TITLE","REMARK","CRYST1","MODEL","ENDMDL"}:
            header.append(ln.rstrip("\n"))
            continue
        if rec != "ATOM":
            continue

        alt = parse_altloc(ln).strip()
        if alt not in keep_altloc:
            # skip altlocs other than blank/'A'
            continue

        resn = parse_resname(ln)
        if resn not in STD_AA:
            # if not a standard AA, drop (for md screening we keep protein only)
            continue

        residues[residue_key(ln)].append(ln.rstrip("\n"))

    out_lines = []
    # Keep header (optional)
    if header:
        out_lines.extend(header)

    dropped = 0
    gly_mutated = 0

    # Write residues in order
    for key in sorted(residues.keys(), key=lambda k: (k[0], k[1], k[2])):
        rlines = residues[key]
        resn = parse_resname(rlines[0])

        atom_names = {parse_atom_name(x) for x in rlines}
        # backbone must exist
        if not {"N","CA","C","O"}.issubset(atom_names):
            dropped += 1
            continue

        # If residue is GLY already, keep as is (but still only ATOM)
        if resn == "GLY":
            out_lines.extend(rlines)
            continue

        # Heuristic for "incomplete sidechain":
        # If residue has ONLY backbone atoms -> it's incomplete already.
        # If residue has some sidechain atoms but is missing others -> also risky.
        # Instead of trying to rebuild, convert to GLY safely.
        # Condition: if any atom name is blank (bad parse) OR sidechain seems truncated
        # We keep it simple: if it doesn't have at least one sidechain heavy atom beyond CB for non-GLY, or if
        # it's known to require many atoms (ARG, LYS, etc.) but has very few.
        sidechain_atoms = atom_names - BACKBONE
        # For non-gly residues, CB normally exists except for PRO sometimes is special; but still has CB.
        # If CB missing or sidechain atom count is suspiciously small, mutate to GLY.
        suspicious = False
        if "CB" not in atom_names:
            suspicious = True
        # Stronger: ARG should have several atoms; if it has fewer than 4 sidechain atoms, it’s likely incomplete.
        if resn == "ARG" and len(sidechain_atoms) < 4:
            suspicious = True
        if resn == "LYS" and len(sidechain_atoms) < 4:
            suspicious = True
        if resn in {"GLU","GLN","MET","PHE","TYR","TRP"} and len(sidechain_atoms) < 3:
            suspicious = True
        if len(sidechain_atoms) == 0:
            suspicious = True

        if suspicious:
            # mutate to GLY: keep only backbone atoms, change resname to GLY
            gly_mutated += 1
            for ln in rlines:
                an = parse_atom_name(ln)
                if an not in BACKBONE:
                    continue
                # replace resname columns 18-20
                ln2 = ln[:17] + "GLY" + ln[20:]
                out_lines.append(ln2)
        else:
            out_lines.extend(rlines)

    out_lines.append("END")

    with open(out_pdb, "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")

    print(f"[sanitize] Wrote: {out_pdb}")
    print(f"[sanitize] Residues mutated to GLY: {gly_mutated}")
    print(f"[sanitize] Residues dropped (missing backbone): {dropped}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("-o","--output", required=True)
    args = ap.parse_args()
    sanitize(args.input, args.output)

if __name__ == "__main__":
    main()
