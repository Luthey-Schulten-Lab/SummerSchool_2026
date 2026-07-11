#!/bin/bash
# Verify setup and submit the overnight DNA Slurm job from the student's home clone.
#
# Run from a Delta login node or OOD Desktop terminal (Tuesday evening):
#   bash ~/SummerSchool_2026/DNA/files/full_cell_simulation/submit_simulation.sh
#
# Or from the instructor copy on bgvl (same effect — submits using ~/SummerSchool_2026):
#   bash /projects/bgvl/SummerSchool_2026/DNA/files/full_cell_simulation/submit_simulation.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${HOME}/SummerSchool_2026"
SIM_DIR="${REPO}/DNA/files/full_cell_simulation"
LOG_DIR="${REPO}/DNA"

bash "${SCRIPT_DIR}/prelaunch_simulation.sh"
mkdir -p "${LOG_DIR}"
sbatch --output="${LOG_DIR}/DNA_tutorial.log" "${SIM_DIR}/launch_simulation.sh"

echo ""
echo "Job submitted. Monitor with:"
echo "  squeue -u ${USER}"
echo "  tail -f ${LOG_DIR}/DNA_tutorial.log"
