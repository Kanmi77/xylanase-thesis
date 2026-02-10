#!/usr/bin/env python3
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_pdb(in_pdb: str, out_pdb: str, ph: float = 7.0) -> None:
    fixer = PDBFixer(filename=in_pdb)

    # Identify missing residues (typically at termini / loops)
    fixer.findMissingResidues()

    # Optional: do NOT build missing internal residues (keeps topology stable for screening)
    # This avoids big loop modeling that can be unpredictable
    if fixer.missingResidues:
        # Remove missing residues that are not at termini (internal gaps)
        chains = list(fixer.topology.chains())
        for chain_idx, chain in enumerate(chains):
            keys_to_delete = []
            for key, residues in fixer.missingResidues.items():
                # key is (chainIndex, residueIndex)
                if key[0] != chain_idx:
                    continue
                # internal missing if not first/last few residues
                # simplest safe rule: keep only terminal missing residues
                # if it’s internal, delete from missingResidues so it won't be built
                res_indices = [r.index for r in chain.residues()]
                if not res_indices:
                    continue
                min_i, max_i = min(res_indices), max(res_indices)
                missing_res_index = key[1]
                if (missing_res_index > min_i + 1) and (missing_res_index < max_i - 1):
                    keys_to_delete.append(key)
            for k in keys_to_delete:
                del fixer.missingResidues[k]

    # Find and rebuild missing heavy atoms (THIS is what fixes missing CD etc.)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Add hydrogens after heavy atoms are complete
    fixer.addMissingHydrogens(pH=ph)

    with open(out_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input PDB")
    ap.add_argument("-o", "--output", required=True, help="Output fixed PDB")
    ap.add_argument("--ph", type=float, default=7.0, help="pH for hydrogenation")
    args = ap.parse_args()
    fix_pdb(args.input, args.output, ph=args.ph)
    print(f"Wrote fixed PDB: {args.output}")

if __name__ == "__main__":
    main()
