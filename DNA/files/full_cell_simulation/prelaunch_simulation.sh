#!/bin/bash
# Verify the student's clone and shared instructor files before submitting.

set -euo pipefail

REPO="${HOME}/SummerSchool_2026"
SIM_DIR="${REPO}/DNA/files/full_cell_simulation"
SHARED_FILES="/projects/bgvl/SummerSchool_2026/DNA/files"

if [[ ! -d "${SIM_DIR}/DNA_SummerSchool_2026/scripts" ]]; then
    echo "ERROR: simulation template not found." >&2
    echo "Clone the repository to ~/SummerSchool_2026 (see DNA/README.md §1):" >&2
    echo "  cd ~ && git clone https://github.com/Luthey-Schulten-Lab/SummerSchool_2026.git --depth 1" >&2
    exit 1
fi

if [[ ! -f "${SIM_DIR}/launch_simulation.sh" ]]; then
    echo "ERROR: ${SIM_DIR}/launch_simulation.sh not found." >&2
    echo "Update your clone: cd ~/SummerSchool_2026 && git pull" >&2
    exit 1
fi

if [[ ! -f "${SHARED_FILES}/DNA_summer2025.sif" ]]; then
    echo "ERROR: ${SHARED_FILES}/DNA_summer2025.sif not found." >&2
    exit 1
fi

if [[ ! -f "${SHARED_FILES}/btree_chromo_gpu/build/apps/btree_chromo" ]]; then
    echo "ERROR: btree_chromo binary not found under ${SHARED_FILES}/btree_chromo_gpu/" >&2
    exit 1
fi

echo "Ready to submit for ${USER}"
echo "  Simulation output: ${SIM_DIR}/DNA_SummerSchool_2026/data/"
echo "  Shared container/binary: ${SHARED_FILES}"
