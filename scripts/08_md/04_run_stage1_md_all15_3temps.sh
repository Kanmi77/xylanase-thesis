#!/usr/bin/env bash
set -euo pipefail

# Run Stage-1 MD for ALL 15 top candidates at 3 temperatures:
# EM -> NVT (100 ps) -> NPT (100 ps) -> MD (1000 ps)
#
# Inputs:
#   results/ranking/top15_final_ranked.csv  (pdb_tag column)
#   docking/proteins/<pdb_tag>.pdb         (repaired PDB from docking prep)
#
# Outputs (per protein, per temp):
#   md/systems/<pdb_tag>/<T>/em.* nvt.* npt.* md.*
#
# Requirements:
#   - GROMACS in conda env: ~/miniconda3/envs/test/bin.AVX2_256/gmx
#   - python scripts:
#       scripts/08_md/00_fix_pdb_for_gromacs.py  (produces input.fixed2.pdb)
#
# Notes:
#   - This script does NOT reuse/copy solv/topol between temperatures; each T builds cleanly to avoid mismatch errors.
#   - Uses amber99sb-ildn + tip3p (as you already used).

BASE="${HOME}/xylanase-thesis"
RANK_CSV="${BASE}/results/ranking/top15_final_ranked.csv"
DOCK_PROT_DIR="${BASE}/docking/proteins"
MD_BASE="${BASE}/md/systems"

GMX="${HOME}/miniconda3/envs/test/bin.AVX2_256/gmx"
export GMX
alias gmx="$GMX"

TEMPS=(333 353 373)

FF="amber99sb-ildn"
WATER="tip3p"
PH="7.0"

# ---- helpers ----
die(){ echo "[FATAL] $*" >&2; exit 1; }
info(){ echo -e "[INFO] $*"; }

need(){
  [[ -e "$1" ]] || die "Missing: $1"
}

write_mdp_em(){
  local out="$1"
  cat > "$out" <<'EOF'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
cutoff-scheme = Verlet
nstlist     = 20
rlist       = 1.1
coulombtype = PME
rcoulomb    = 1.1
rvdw        = 1.1
pbc         = xyz
EOF
}

write_mdp_nvt(){
  local out="$1"
  local tempK="$2"
  cat > "$out" <<EOF
integrator  = md
dt          = 0.002
nsteps      = 50000        ; 100 ps
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500

cutoff-scheme = Verlet
nstlist     = 50
rlist       = 1.1
coulombtype = PME
rcoulomb    = 1.1
rvdw        = 1.1

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = ${tempK} ${tempK}

pcoupl      = no
constraints = h-bonds
constraint_algorithm = lincs

pbc         = xyz
gen_vel     = yes
gen_temp    = ${tempK}
gen_seed    = -1
EOF
}

write_mdp_npt(){
  local out="$1"
  local tempK="$2"
  cat > "$out" <<EOF
integrator  = md
dt          = 0.002
nsteps      = 50000        ; 100 ps
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500

cutoff-scheme = Verlet
nstlist     = 50
rlist       = 1.1
coulombtype = PME
rcoulomb    = 1.1
rvdw        = 1.1

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = ${tempK} ${tempK}

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

constraints = h-bonds
constraint_algorithm = lincs
pbc         = xyz
gen_vel     = no
EOF
}

write_mdp_md(){
  local out="$1"
  local tempK="$2"
  cat > "$out" <<EOF
integrator  = md
dt          = 0.002
nsteps      = 500000       ; 1000 ps = 1 ns
nstxout-compressed = 1000
nstenergy   = 1000
nstlog      = 1000

cutoff-scheme = Verlet
nstlist     = 50
rlist       = 1.1
coulombtype = PME
rcoulomb    = 1.1
rvdw        = 1.1

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1
ref_t       = ${tempK} ${tempK}

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

constraints = h-bonds
constraint_algorithm = lincs
pbc         = xyz
gen_vel     = no
EOF
}

prep_one_temp(){
  local pdb_tag="$1"
  local tempK="$2"

  local prot_dir="${MD_BASE}/${pdb_tag}/${tempK}K"
  mkdir -p "$prot_dir"
  info "Protein: ${pdb_tag} | Temp: ${tempK}K -> ${prot_dir}"

  # 0) Input PDB
  local src_pdb="${DOCK_PROT_DIR}/${pdb_tag}.pdb"
  need "$src_pdb"
  cp -f "$src_pdb" "${prot_dir}/input.pdb"

  # 1) Fix PDB (ensures ARG CD etc)
  python "${BASE}/scripts/08_md/00_fix_pdb_for_gromacs.py" \
    -i "${prot_dir}/input.pdb" \
    -o "${prot_dir}/input.fixed2.pdb" \
    --ph "${PH}"

  # 2) Topology
  ( cd "$prot_dir" && \
    $GMX pdb2gmx -f input.fixed2.pdb -o processed.gro -p topol.top -i posre.itp -ff "$FF" -water "$WATER" -ignh \
  )

  # 3) Box + solvate
  ( cd "$prot_dir" && \
    $GMX editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic && \
    $GMX solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top \
  )

  # 4) Ions
  ( cd "$prot_dir" && \
    write_mdp_em ions.mdp && \
    $GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 && \
    printf "SOL\n" | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 \
  )

  # 5) EM
  ( cd "$prot_dir" && \
    write_mdp_em em.mdp && \
    $GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1 && \
    $GMX mdrun -deffnm em \
  )

  # 6) NVT
  ( cd "$prot_dir" && \
    write_mdp_nvt nvt.mdp "$tempK" && \
    $GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1 && \
    $GMX mdrun -deffnm nvt \
  )

  # 7) NPT
  ( cd "$prot_dir" && \
    write_mdp_npt npt.mdp "$tempK" && \
    $GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 && \
    $GMX mdrun -deffnm npt \
  )

  # 8) Production MD
  ( cd "$prot_dir" && \
    write_mdp_md md.mdp "$tempK" && \
    $GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 1 && \
    $GMX mdrun -deffnm md \
  )

  info "DONE: ${pdb_tag} ${tempK}K"
}

main(){
  need "$RANK_CSV"
  need "$GMX"
  need "${BASE}/scripts/08_md/00_fix_pdb_for_gromacs.py"

  mkdir -p "$MD_BASE"

  info "Reading: $RANK_CSV"
  # Extract pdb_tag column from CSV robustly (expects header with pdb_tag)
  mapfile -t TAGS < <(python - <<'PY'
import csv, os
base=os.path.expanduser("~/xylanase-thesis")
path=os.path.join(base,"results/ranking/top15_final_ranked.csv")
with open(path,newline="") as f:
    r=csv.DictReader(f)
    if "pdb_tag" not in r.fieldnames:
        raise SystemExit("Missing column pdb_tag in top15_final_ranked.csv")
    tags=[row["pdb_tag"].strip() for row in r if row.get("pdb_tag")]
print("\n".join(tags))
PY
)

  info "Found ${#TAGS[@]} proteins"
  [[ "${#TAGS[@]}" -eq 15 ]] || info "Note: expected 15, got ${#TAGS[@]} (still proceeding)."

  for tag in "${TAGS[@]}"; do
    for T in "${TEMPS[@]}"; do
      prep_one_temp "$tag" "$T"
    done
  done

  info "ALL DONE: Stage-1 MD for all tags + all temperatures."
}

main "$@"
