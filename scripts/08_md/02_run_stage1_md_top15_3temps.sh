#!/usr/bin/env bash
set -euo pipefail

# Path Configuration
BASE="$HOME/xylanase-thesis"
INP="$BASE/md/inputs"
SYS="$BASE/md/systems"
MDP="$BASE/md/mdp"
LOG="$BASE/md/logs"

mkdir -p "$SYS" "$LOG"

# Stage 1 temperatures (333K=60°C, 353K=80°C, 373K=100°C)
TEMPS=("333" "353" "373")

# Force field settings (AMBER99SB-ILDN is excellent for side-chain torsion accuracy)
FF="amber99sb-ildn"
WATER="tip3p"

# Check for input files
PDBS=("$INP/1O8S_A_Repair.pdb")
if [ ${#PDBS[@]} -eq 0 ]; then
  echo "ERROR: No PDBs found in $INP"
  exit 1
fi

echo "[INFO] Found ${#PDBS[@]} input PDBs"
echo "[INFO] Temperatures: ${TEMPS[*]} K"
echo "[INFO] Force field: $FF | Water: $WATER"

for pdb in "${PDBS[@]}"; do
  tag=$(basename "$pdb" .pdb)
  echo ""
  echo "============================================================"
  echo "[INFO] Processing Protein: $tag"
  echo "============================================================"

  for T in "${TEMPS[@]}"; do
    wd="$SYS/$tag/${T}K"
    mkdir -p "$wd"
    echo "[INFO] Current Target: ${T}K -> $wd"

    # Step 0: Transfer and Sanitize PDB
    cp -f "$pdb" "$wd/input.pdb"

    # Fix PDB: This is the "Nuclear Option" to remove Cd, fix missing atoms, 
    # and set protonation states for pH 7.0 (Master Thesis standard)
    echo "[INFO] Step 0: Fixing PDB for GROMACS..."
    python "$BASE/scripts/08_md/00_fix_pdb_for_gromacs.py" \
      -i "$wd/input.pdb" \
      -o "$wd/input.fixed.pdb" \
      --ph 7.0 2>&1 | tee "$LOG/${tag}_${T}K_pdbfixer.log"

    # Step 1: Generate Topology
    # Using 'input.fixed.pdb' ensures pdb2gmx won't crash on missing atoms
    echo "[INFO] Step 1: Generating Topology (pdb2gmx)..."
    gmx pdb2gmx -f "$wd/input.fixed.pdb" -o "$wd/processed.gro" -p "$wd/topol.top" -i "$wd/posre.itp" \
      -ff "$FF" -water "$WATER" -ignh 2>&1 | tee "$LOG/${tag}_${T}K_pdb2gmx.log"

    # Step 2: Create Simulation Box
    echo "[INFO] Step 2: Defining Box (editconf)..."
    gmx editconf -f "$wd/processed.gro" -o "$wd/newbox.gro" -c -d 1.0 -bt cubic \
      2>&1 | tee "$LOG/${tag}_${T}K_editconf.log"

    # Step 3: Solvate with Water
    echo "[INFO] Step 3: Solvating (solvate)..."
    gmx solvate -cp "$wd/newbox.gro" -cs spc216.gro -o "$wd/solv.gro" -p "$wd/topol.top" \
      2>&1 | tee "$LOG/${tag}_${T}K_solvate.log"

    # Step 4: Neutralize and Add Ions (0.15 M NaCl is physiological/standard)
    echo "[INFO] Step 4: Adding Ions (genion)..."
    gmx grompp -f "$MDP/minim.mdp" -c "$wd/solv.gro" -p "$wd/topol.top" -o "$wd/ions.tpr" \
      -maxwarn 2 2>&1 | tee "$LOG/${tag}_${T}K_grompp_ions.log"
    
    echo "SOL" | gmx genion -s "$wd/ions.tpr" -o "$wd/solv_ions.gro" -p "$wd/topol.top" \
      -pname NA -nname CL -neutral -conc 0.15 \
      2>&1 | tee "$LOG/${tag}_${T}K_genion.log"

    # Step 5: Energy Minimization (Relieve steric clashes)
    echo "[INFO] Step 5: Energy Minimization..."
    gmx grompp -f "$MDP/minim.mdp" -c "$wd/solv_ions.gro" -p "$wd/topol.top" -o "$wd/em.tpr" \
      -maxwarn 2 2>&1 | tee "$LOG/${tag}_${T}K_grompp_em.log"
    gmx mdrun -deffnm "$wd/em" 2>&1 | tee "$LOG/${tag}_${T}K_mdrun_em.log"

    # Step 6: NVT Equilibration (Heat up to target temperature)
    echo "[INFO] Step 6: NVT Equilibration at ${T}K..."
    sed "s/TEMP/${T}/g" "$MDP/nvt.mdp" > "$wd/nvt_${T}.mdp"
    gmx grompp -f "$wd/nvt_${T}.mdp" -c "$wd/em.gro" -r "$wd/em.gro" -p "$wd/topol.top" -o "$wd/nvt.tpr" \
      -maxwarn 2 2>&1 | tee "$LOG/${tag}_${T}K_grompp_nvt.log"
    gmx mdrun -deffnm "$wd/nvt" 2>&1 | tee "$LOG/${tag}_${T}K_mdrun_nvt.log"

    # Step 7: NPT Equilibration (Stabilize Pressure)
    echo "[INFO] Step 7: NPT Equilibration..."
    sed "s/TEMP/${T}/g" "$MDP/npt.mdp" > "$wd/npt_${T}.mdp"
    gmx grompp -f "$wd/npt_${T}.mdp" -c "$wd/nvt.gro" -r "$wd/nvt.gro" -t "$wd/nvt.cpt" -p "$wd/topol.top" -o "$wd/npt.tpr" \
      -maxwarn 2 2>&1 | tee "$LOG/${tag}_${T}K_grompp_npt.log"
    gmx mdrun -deffnm "$wd/npt" 2>&1 | tee "$LOG/${tag}_${T}K_mdrun_npt.log"

    # Step 8: Production MD (The actual data collection)
    echo "[INFO] Step 8: Starting Production MD..."
    sed "s/TEMP/${T}/g" "$MDP/md.mdp" > "$wd/md_${T}.mdp"
    gmx grompp -f "$wd/md_${T}.mdp" -c "$wd/npt.gro" -t "$wd/npt.cpt" -p "$wd/topol.top" -o "$wd/md.tpr" \
      -maxwarn 2 2>&1 | tee "$LOG/${tag}_${T}K_grompp_md.log"
    
    # Use -cpi to resume if the server goes down
    gmx mdrun -deffnm "$wd/md" -cpi "$wd/md.cpt" 2>&1 | tee "$LOG/${tag}_${T}K_mdrun_md.log"

    echo "[SUCCESS] Finished Pipeline for ${tag} at ${T}K"
  done
done

echo ""
echo "[FINISH] All proteins simulated at all temperatures."
