#!/bin/bash
#SBATCH --account=bgvl-delta-gpu
#SBATCH --job-name=DNA_test
#SBATCH --partition=gpuA100x4
#SBATCH --time=02:00:00
#SBATCH --mem=32g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
#SBATCH --nodes=1

# 2-minute biological-time smoke test (timesteps 0–1).
#
#   sbatch --output=~/SummerSchool_2026/DNA/DNA_test.log \
#     ~/SummerSchool_2026/DNA/files/full_cell_simulation/launch_test.sh

set -euo pipefail

SIM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_FILES="/projects/bgvl/SummerSchool_2026/DNA/files"
FILES_ROOT="${SHARED_FILES}"

if [[ -f "${SIM_DIR}/../DNA_summer2025.sif" ]]; then
  FILES_ROOT="$(cd "${SIM_DIR}/.." && pwd)"
fi

SIF="${FILES_ROOT}/DNA_summer2025.sif"
BTREE_BIN="${FILES_ROOT}/btree_chromo_gpu/build/apps/btree_chromo"
TEMPLATE_DIR="${SIM_DIR}/DNA_SummerSchool_2026"
SIM_ROOT="${1:-${SIM_DIR}/DNA_SummerSchool_2026_test}"

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: ${SIF} not found." >&2
  exit 1
fi

if [[ ! -f "${BTREE_BIN}" ]]; then
  echo "ERROR: ${BTREE_BIN} not found. Run: bash ${SHARED_FILES}/build_btree_chromo.sh" >&2
  exit 1
fi

mkdir -p "${SIM_ROOT}"
echo "Staging simulation template into ${SIM_ROOT}"
cp -rp "${TEMPLATE_DIR}/." "${SIM_ROOT}/"

apptainer run \
  --nv \
  --writable-tmpfs \
  --no-home \
  --containall \
  --bind "${SIM_ROOT}:/mnt" \
  --bind "${FILES_ROOT}:/ps:ro" \
  "${SIF}" /bin/bash -c "cd /mnt/scripts && bash run_sc_chain_generation.sh && python3 run_btree_chromo.py 34 summerschool 0 1"
