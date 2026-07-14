#!/bin/bash
#SBATCH --account=bgvl-delta-gpu
#SBATCH --job-name=DNA_tutorial
#SBATCH --partition=gpuA100x4
#SBATCH --time=15:00:00
#SBATCH --mem=32g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
#SBATCH --nodes=1

# Full-cell simulation using shared DNA_summer2025.sif and btree_chromo on bgvl.
# Simulation reads/writes in-place under ~/SummerSchool_2026/.../DNA_SummerSchool_2026/.
#
# Submit from a Delta login node or OOD Desktop terminal:
#   bash ~/SummerSchool_2026/DNA/files/full_cell_simulation/submit_simulation.sh

set -euo pipefail

# Slurm copies this script to /var/spool/slurmd/job*/ — do not resolve paths from BASH_SOURCE.
REPO="${HOME}/SummerSchool_2026"
SIM_DIR="${REPO}/DNA/files/full_cell_simulation"
SIM_ROOT="${SIM_DIR}/DNA_SummerSchool_2026"
SHARED_FILES="/projects/bgvl/SummerSchool_2026/DNA/files"
FILES_ROOT="${SHARED_FILES}"

if [[ -f "${SIM_DIR}/../DNA_summer2025.sif" ]]; then
  FILES_ROOT="$(cd "${SIM_DIR}/.." && pwd)"
fi

SIF="${FILES_ROOT}/DNA_summer2025.sif"
BTREE_BIN="${FILES_ROOT}/btree_chromo_gpu/build/apps/btree_chromo"

if [[ ! -d "${SIM_ROOT}/scripts" ]]; then
  echo "ERROR: ${SIM_ROOT}/scripts not found." >&2
  echo "Clone the repository to ~/SummerSchool_2026 (see DNA/README.md §1)." >&2
  exit 1
fi

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: ${SIF} not found." >&2
  exit 1
fi

if [[ ! -f "${BTREE_BIN}" ]]; then
  echo "ERROR: ${BTREE_BIN} not found." >&2
  echo "Run: bash ${SHARED_FILES}/build_btree_chromo.sh" >&2
  exit 1
fi

# Per-student RNG seed (stable for a given username; override with DNA_SIM_SEED if needed).
if [[ -n "${DNA_SIM_SEED:-}" ]]; then
  SIM_SEED="${DNA_SIM_SEED}"
else
  SIM_SEED="$(
    python3 -c "import hashlib; u='${USER:-unknown}'; print(int(hashlib.sha256(u.encode()).hexdigest()[:7], 16) % 90 + 10)"
  )"
fi
echo "Simulation RNG seed for ${USER}: ${SIM_SEED}"

# Open OnDemand Desktop nests Apptainer with binds like $OOD_CWD:/cwd and
# $TMPDIR:/TMPDIR. If those appear in APPTAINER_BINDPATH / SINGULARITY_BINDPATH
# when the job is submitted from an OOD terminal, the compute-node Apptainer
# tries to mount host paths /cwd and /TMPDIR (which do not exist) and fails.
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH || true
unset APPTAINER_TMPDIR SINGULARITY_TMPDIR || true

# Mount DNA/files at /ps so scripts see /ps/btree_chromo_gpu/...
apptainer run \
  --nv \
  --writable-tmpfs \
  --no-home \
  --containall \
  --pwd /mnt/scripts \
  --bind "${SIM_ROOT}:/mnt" \
  --bind "${FILES_ROOT}:/ps:ro" \
  "${SIF}" /bin/bash -c "export SIM_SEED=${SIM_SEED}; bash run_sc_chain_generation.sh && python3 run_btree_chromo.py ${SIM_SEED} summerschool 0 90"
