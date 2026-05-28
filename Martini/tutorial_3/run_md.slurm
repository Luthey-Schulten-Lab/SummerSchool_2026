#!/bin/bash
#SBATCH --job-name=bentopy_md
#SBATCH --account=bgvl-delta-gpu
#SBATCH --partition=gpuH200x8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --output=md_%j.out
#SBATCH --error=md_%j.err

set -euo pipefail

module purge
module load gromacs/2026.1-cuda

cd "$SLURM_SUBMIT_DIR"
echo "Host:        $(hostname)"
echo "Submit dir:  $SLURM_SUBMIT_DIR"
echo "GPUs:        ${CUDA_VISIBLE_DEVICES:-unknown}"
gmx_mpi --version | head -5
echo "================================================================"

NCPU="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="$NCPU"

MDRUN="srun --ntasks=1 --cpus-per-task=$NCPU gmx_mpi mdrun -nb gpu -pin on -ntomp $NCPU -v"

# ---------------- Energy minimization ----------------
if [[ ! -s em.gro ]]; then
    echo "[$(date)] >>> Energy minimization"
    gmx_mpi grompp -f mdp_files/em.mdp -c solvated_system.gro -p topol.top \
        -o em.tpr -maxwarn 1
    $MDRUN -deffnm em
else
    echo "[$(date)] em.gro exists, skipping EM"
fi

# ---------------- Index file (Protein / Lipid / Solvent) ----------------
if [[ ! -s index.ndx ]]; then
    echo "[$(date)] >>> Building index.ndx"
    gmx_mpi make_ndx -f em.gro -o index.ndx <<EOF
name 13 Lipid
r W | r ION
name 16 Solvent
q
EOF
else
    echo "[$(date)] index.ndx exists, skipping make_ndx"
fi

# ---------------- Equilibration (100 ps) ----------------
if [[ ! -s eq.gro ]]; then
    echo "[$(date)] >>> Equilibration"
    gmx_mpi grompp -f mdp_files/eq.mdp -c em.gro -p topol.top -n index.ndx \
        -o eq.tpr -maxwarn 1
    $MDRUN -deffnm eq
else
    echo "[$(date)] eq.gro exists, skipping equilibration"
fi

# ---------------- Production (1 us) with checkpoint resume ----------------
if [[ ! -s md.tpr ]]; then
    echo "[$(date)] >>> grompp production"
    gmx_mpi grompp -f mdp_files/md.mdp -c eq.gro -p topol.top -n index.ndx \
        -o md.tpr -maxwarn 1
fi

# -maxh 47.5 leaves a small buffer before the 48 h time limit so GROMACS
# writes a clean checkpoint and exits gracefully.
if [[ -s md.cpt ]]; then
    echo "[$(date)] >>> Resuming production from md.cpt"
    $MDRUN -deffnm md -cpi md.cpt -append -maxh 47.5
else
    echo "[$(date)] >>> Starting production from scratch"
    $MDRUN -deffnm md -maxh 47.5
fi

echo "[$(date)] >>> Job script complete"
