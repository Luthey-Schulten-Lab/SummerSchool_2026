#!/usr/bin/env bash
# Copy 4DWCM source + launch script into /projects/bgvl/$USER/4dwcm_run/
set -euo pipefail

GATEWAY="/projects/bgvl/containers/4DWCM_Gateway"
RUN_DIR="/projects/bgvl/${USER}/4dwcm_run"

if [[ ! -f "${GATEWAY}/Optimize_4DWCM_Minimal_Cell/Whole_Cell_Minimal_Cell.py" ]]; then
  echo "ERROR: source not found under ${GATEWAY}/Optimize_4DWCM_Minimal_Cell"
  exit 1
fi
if [[ ! -f "${GATEWAY}/4DWCM_ssh/launch_4dwcm_7200.sh" ]]; then
  echo "ERROR: scripts not found under ${GATEWAY}/4DWCM_ssh"
  exit 1
fi

mkdir -p "${RUN_DIR}/logs"

echo "=== Copy 4DWCM into your folder ==="
echo "Your folder: ${RUN_DIR}"
echo ""

echo "Copying Optimize_4DWCM_Minimal_Cell ..."
cp -r "${GATEWAY}/Optimize_4DWCM_Minimal_Cell" "${RUN_DIR}/"

echo "Copying launch scripts ..."
cp -f "${GATEWAY}/4DWCM_ssh/launch_4dwcm_7200.sh" \
      "${GATEWAY}/4DWCM_ssh/launch_4dwcm_restart.sh" \
      "${RUN_DIR}/"
chmod +x "${RUN_DIR}/launch_4dwcm_7200.sh" "${RUN_DIR}/launch_4dwcm_restart.sh"

echo ""
echo "Done."
echo ""
echo "Submit job:"
echo "  cd ${RUN_DIR}"
echo "  sbatch launch_4dwcm_7200.sh"
echo ""
echo "Outputs will be written under YOUR folder (not 4DWCM_Gateway):"
echo "  Slurm log  : ${RUN_DIR}/logs/4dwcm_7200-<JOBID>.out"
echo "  Science data: ${RUN_DIR}/Optimize_4DWCM_Minimal_Cell/Data/4dwcm_7200/"
