#!/usr/bin/env python3
"""
07_docking/01_prepare_docking_inputs.py

Prepare docking inputs for top-15 candidates:
- takes FoldX-repaired PDBs (foldx/repair/<pdb_tag>.pdb)
- cleans receptor PDB (protein-only, remove waters/ligands/metals, altloc handling)
- fixes PDB element column (needed by meeko/prody/rdkit)
- generates receptor PDBQT via meeko (mk_prepare_receptor.py) using:
    --allow_bad_res --default_altloc A -o <basename> --write_pdbqt
- writes Vina config files with a simple centroid-based box (CA centroid fallback to all atoms)
- writes a manifest CSV for downstream Vina batch runs

This version matches your Meeko behavior: "prepared but no files written" unless --write_pdbqt is provided.
"""

import os
import sys
import subprocess
from typing import Tuple, List, Dict

import numpy as np
import pandas as pd


BASE = os.path.expanduser("~/xylanase-thesis")

TOP15_CSV = os.path.join(BASE, "results/ranking/top15_candidates.csv")
REPAIR_DIR = os.path.join(BASE, "foldx/repair")

DOCK_DIR = os.path.join(BASE, "docking")
PROT_DIR = os.path.join(DOCK_DIR, "proteins")
CONF_DIR = os.path.join(DOCK_DIR, "configs")
LOG_DIR = os.path.join(DOCK_DIR, "logs")
MANIFEST_DIR = os.path.join(DOCK_DIR, "manifests")

OUT_MANIFEST = os.path.join(MANIFEST_DIR, "dock_manifest_top15.csv")

# These are scripts you created earlier (must exist)
CLEANER = os.path.join(BASE, "scripts/07_docking/00_clean_receptor_pdb.py")
FIXER = os.path.join(BASE, "scripts/07_docking/00_fix_pdb_elements.py")

# Vina box sizes: keep conservative; GH10 larger than GH11
BOX_SIZES = {
    "GH10": (28.0, 28.0, 28.0),
    "GH11": (22.0, 22.0, 22.0),
}
DEFAULT_BOX = (24.0, 24.0, 24.0)


def ensure_dirs() -> None:
    for d in [DOCK_DIR, PROT_DIR, CONF_DIR, LOG_DIR, MANIFEST_DIR]:
        os.makedirs(d, exist_ok=True)


def run(cmd: List[str], log_path: str = "") -> str:
    """Run a subprocess, capture stdout+stderr, optionally log output, raise on failure."""
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = r.stdout
    if log_path:
        with open(log_path, "w", encoding="utf-8") as f:
            f.write(out)
    if r.returncode != 0:
        raise RuntimeError(f"Command failed rc={r.returncode}: {' '.join(cmd)}\n{out}")
    return out


def which_or_fail(exe: str) -> None:
    """Fail early if executable is not in PATH."""
    if subprocess.run(["bash", "-lc", f"command -v {exe} >/dev/null 2>&1"]).returncode != 0:
        raise RuntimeError(f"Required executable not found in PATH: {exe}")


def centroid_from_pdb(pdb_path: str) -> Tuple[float, float, float]:
    """
    Compute centroid of CA atoms if present, else fallback to all ATOM/HETATM coordinates.
    """
    ca = []
    atoms = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            # coords in PDB fixed columns
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except Exception:
                continue
            atoms.append((x, y, z))
            atom_name = line[12:16].strip()
            if atom_name == "CA":
                ca.append((x, y, z))
    coords = ca if ca else atoms
    if not coords:
        raise ValueError(f"No coordinates parsed from {pdb_path}")
    arr = np.array(coords, dtype=float)
    c = arr.mean(axis=0)
    return float(c[0]), float(c[1]), float(c[2])


def make_vina_conf(conf_path: str, receptor_pdbqt: str,
                  center: Tuple[float, float, float],
                  size: Tuple[float, float, float]) -> None:
    cx, cy, cz = center
    sx, sy, sz = size
    with open(conf_path, "w", encoding="utf-8") as f:
        f.write(f"receptor = {receptor_pdbqt}\n")
        f.write(f"center_x = {cx:.3f}\n")
        f.write(f"center_y = {cy:.3f}\n")
        f.write(f"center_z = {cz:.3f}\n")
        f.write(f"size_x = {sx:.1f}\n")
        f.write(f"size_y = {sy:.1f}\n")
        f.write(f"size_z = {sz:.1f}\n")
        f.write("exhaustiveness = 12\n")
        f.write("num_modes = 9\n")


def prepare_one(pdb_tag: str, gh_family: str, meta: Dict[str, str]) -> Dict[str, object]:
    """
    Prepare one receptor:
      - copy repaired PDB
      - clean -> fixed
      - meeko receptor pdbqt using basename output + --write_pdbqt
      - compute centroid and write vina config
    """
    repaired_in = os.path.join(REPAIR_DIR, f"{pdb_tag}.pdb")
    if not os.path.exists(repaired_in):
        raise FileNotFoundError(f"Missing repaired PDB: {repaired_in}")

    # Keep a local copy in docking/proteins
    repaired_copy = os.path.join(PROT_DIR, f"{pdb_tag}.pdb")
    if not os.path.exists(repaired_copy):
        with open(repaired_in, "rb") as src, open(repaired_copy, "wb") as dst:
            dst.write(src.read())

    # Clean receptor (protein only, remove hetero atoms & metals)
    protein_clean = os.path.join(PROT_DIR, f"{pdb_tag}.protein_clean.pdb")
    clean_log = os.path.join(LOG_DIR, f"{pdb_tag}.clean.log")
    if not os.path.exists(protein_clean):
        run(["python", CLEANER, repaired_copy, protein_clean], log_path=clean_log)

    # Fix element column
    fixed_pdb = os.path.join(PROT_DIR, f"{pdb_tag}.protein_clean.fixed.pdb")
    fix_log = os.path.join(LOG_DIR, f"{pdb_tag}.fix_elements.log")
    if not os.path.exists(fixed_pdb):
        run(["python", FIXER, protein_clean, fixed_pdb], log_path=fix_log)

    # Meeko receptor PDBQT:
    # IMPORTANT: your Meeko needs -o BASENAME plus --write_pdbqt
    out_base = os.path.join(PROT_DIR, pdb_tag)  # basename, no extension
    receptor_pdbqt = out_base + ".pdbqt"
    prep_log = os.path.join(LOG_DIR, f"{pdb_tag}.receptor_prep.log")

    if not os.path.exists(receptor_pdbqt):
        cmd = [
            "mk_prepare_receptor.py",
            "--allow_bad_res",
            "--default_altloc", "A",
            "-i", fixed_pdb,
            "-o", out_base,
            "--write_pdbqt",
        ]
        run(cmd, log_path=prep_log)

    # Validate output exists
    if not os.path.exists(receptor_pdbqt) or os.path.getsize(receptor_pdbqt) == 0:
        raise RuntimeError(
            f"Receptor PDBQT not produced for {pdb_tag}: {receptor_pdbqt}\n"
            f"Check: {prep_log}"
        )

    # Vina config
    center = centroid_from_pdb(fixed_pdb)
    size = BOX_SIZES.get(gh_family, DEFAULT_BOX)
    conf_path = os.path.join(CONF_DIR, f"{pdb_tag}.vina.conf")
    make_vina_conf(conf_path, receptor_pdbqt, center, size)

    out = {
        "pdb_tag": pdb_tag,
        "gh_family": gh_family,
        "repaired_pdb": repaired_copy,
        "protein_clean_pdb": protein_clean,
        "fixed_pdb": fixed_pdb,
        "receptor_pdbqt": receptor_pdbqt,
        "vina_conf": conf_path,
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": size[0],
        "size_y": size[1],
        "size_z": size[2],
    }
    # add metadata columns if present
    out.update(meta)
    return out


def main() -> int:
    ensure_dirs()

    for p in [TOP15_CSV, REPAIR_DIR, CLEANER, FIXER]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required path: {p}")

    which_or_fail("mk_prepare_receptor.py")

    df = pd.read_csv(TOP15_CSV, dtype=str).fillna("")
    if "pdb_tag" not in df.columns:
        raise ValueError(f"Expected column pdb_tag in {TOP15_CSV}. Found: {df.columns.tolist()}")
    if "gh_family" not in df.columns:
        raise ValueError(f"Expected column gh_family in {TOP15_CSV}. Found: {df.columns.tolist()}")

    rows = []
    for _, r in df.iterrows():
        pdb_tag = r["pdb_tag"].strip()
        gh = r["gh_family"].strip() if r["gh_family"].strip() else "unknown"

        # keep useful metadata for downstream analysis
        meta = {
            "pdb_id": r.get("pdb_id", ""),
            "chosen_chain": r.get("chosen_chain", ""),
            "organism_type": r.get("organism_type", ""),
            "organism": r.get("organism", ""),
            "uniprot_accession": r.get("uniprot_accession", ""),
        }

        print(f"[prepare] {pdb_tag} ({gh})")
        out = prepare_one(pdb_tag, gh, meta)
        rows.append(out)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_MANIFEST, index=False)
    print(f"Wrote: {OUT_MANIFEST} rows={len(out_df)}")

    # helpful counts
    pdbqt_count = len([x for x in os.listdir(PROT_DIR) if x.endswith(".pdbqt")])
    conf_count = len([x for x in os.listdir(CONF_DIR) if x.endswith(".vina.conf")])
    print(f"PDBQT count: {pdbqt_count}  |  Vina conf count: {conf_count}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
